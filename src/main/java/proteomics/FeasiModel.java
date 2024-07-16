package proteomics;

import ProteomicsLibrary.*;
import org.apache.commons.math3.util.Pair;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.PTM.InferPTM;
import proteomics.Types.*;
import proteomics.Index.BuildIndex;
import proteomics.Parameter.Parameter;
import proteomics.Spectrum.DatasetReader;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;


import java.io.*;
import java.net.InetAddress;
import java.net.UnknownHostException;
import java.sql.*;
import java.text.DecimalFormat;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.locks.ReentrantLock;

public class PIPI {
    private static final Logger logger = LoggerFactory.getLogger(PIPI.class);
    public static final boolean nTermSpecific = true;
    public static final boolean cTermSpecific = true;
    public static final double MIN_PEAK_SUM_INFER_AA = 0.0;
    public static HashSet<Integer> debugScanNum = new HashSet<>(Arrays.asList(123184));//129543, 111179, 109395
    public static void main(String[] args) {
        long startTime = System.nanoTime();
        // Set parameters
        String parameterPath = args[0].trim();
        logger.info("Integer Linear Programming Model For PTM Characterization.");

        String dbName = null;
        String hostName = "unknown-host";
        try {
            hostName = InetAddress.getLocalHost().getHostName();
        } catch (UnknownHostException ex) {
            logger.warn("Cannot get the computer's name.");
        }

        try {
            logger.info("parameter: {}.", parameterPath);
            dbName = String.format(Locale.US, "PIPI.%s.temp.db", hostName);
            new PIPI(parameterPath, dbName, hostName);
        } catch (Exception ex) {
            ex.printStackTrace();
            logger.error(ex.toString());
        } finally {
            if (dbName != null) {
                (new File(dbName)).delete();
                (new File(dbName + "-wal")).delete();
                (new File(dbName + "-shm")).delete();
            }
        }

        double totalHour = (double) (System.nanoTime() - startTime) * 1e-9 / 3600;
        logger.info("Running time: {} hours.", totalHour);
        logger.info("Done!");
    }

    private PIPI(String parameterPath, String dbName, String hostName) throws Exception {
        // Get the parameter map
        Parameter parameter = new Parameter(parameterPath);
        Map<String, String> parameterMap = parameter.returnParameterMap();
        double ms2Tolerance = Double.valueOf(parameterMap.get("ms2_tolerance"));
        double ms1Tolerance = Double.valueOf(parameterMap.get("ms1_tolerance"));
        int ms1ToleranceUnit = Integer.valueOf(parameterMap.get("ms1_tolerance_unit"));
        double minClear = 112.5;
        double maxClear = 121.5;
        String truthPath = parameterMap.get("truth_path");
        String spectraPath = parameterMap.get("spectra_path");
        String outputDir = parameterMap.get("output_dir");

        Set<Integer> msLevelSet = new HashSet<>();
        msLevelSet.add(2);
        logger.info("Loading parameters...");
        BuildIndex buildIndex = new BuildIndex(parameterMap);
        MassTool massTool = buildIndex.returnMassTool();
        InferPTM inferPTM = buildIndex.getInferPTM();

        logger.info("Reading spectra...");
        File spectraFile = new File(spectraPath);
        DatasetReader datasetReader;
        JMzReader[] spectraParserArray;
        String sqlPath = "jdbc:sqlite:" + dbName;
        Class.forName("org.sqlite.JDBC").newInstance();
        Map<Integer, String> fileIdNameMap = new HashMap<>();
        Map<String, Integer> fileNameIdMap = new HashMap<>();
        if ((!spectraFile.exists())) {
            throw new FileNotFoundException("The spectra file not found.");
        }

        if ( ! spectraFile.isDirectory()) {
            spectraParserArray = new JMzReader[1];
            JMzReader spectraParser;
            String ext = spectraPath.substring(spectraPath.lastIndexOf(".")+1);
            spectraParser = new MgfFile(spectraFile);
            spectraParserArray[0] = spectraParser;
            fileIdNameMap.put(0, spectraPath.substring(spectraPath.lastIndexOf("/")+1).split("\\.")[0].replaceAll("\\.","_"));
            fileNameIdMap.put(spectraPath.substring(spectraPath.lastIndexOf("/")+1).split("\\.")[0].replaceAll("\\.","_"), 0);
            datasetReader = new DatasetReader(spectraParserArray, ms1Tolerance, ms1ToleranceUnit, massTool, ext, msLevelSet, sqlPath, fileIdNameMap);
        } else {
            String[] fileList = spectraFile.list(new FilenameFilter() {
                @Override
                public boolean accept(File dir, String name) {
                    return name.endsWith(".mgf");
                }
            });
            spectraParserArray = new JMzReader[fileList.length];
            for (int i = 0; i < fileList.length; i++){
                spectraParserArray[i] = new MgfFile(new File(spectraPath + fileList[i]));
                fileIdNameMap.put(i, fileList[i].split("\\.")[0].replaceAll("\\.","_"));
                fileNameIdMap.put(fileList[i].split("\\.")[0].replaceAll("\\.","_"), i);
            }

            String ext = fileList[0].substring(fileList[0].lastIndexOf(".")+1);
            datasetReader = new DatasetReader(spectraParserArray, ms1Tolerance, ms1ToleranceUnit, massTool, ext, msLevelSet, sqlPath, fileIdNameMap);
        }
        BufferedReader parameterReader = new BufferedReader(new FileReader(truthPath));
        Map<Integer, List<String>> scanNum_truth_map = new HashMap<>();
        String line;
        while ((line = parameterReader.readLine()) != null) {
            if (line.isEmpty()) continue;
            if (line.contains("scanNo")) continue;
            line = line.trim();
            String[] splitRes = line.split(",");
            scanNum_truth_map.put(Integer.valueOf(splitRes[0]), Arrays.asList(splitRes).subList(1, splitRes.length));
        }

        SpecProcessor specProcessor = new SpecProcessor(massTool);

        logger.info("Starts...");
        int threadNum = Integer.valueOf(parameterMap.get("thread_num"));
        if (threadNum == 0) {
            threadNum = 3 + Runtime.getRuntime().availableProcessors();
        }
        if (java.lang.management.ManagementFactory.getRuntimeMXBean().getInputArguments().toString().indexOf("jdwp") >= 0){
            threadNum = 1;
        }

        //////////==================================================
        ExecutorService threadPool = Executors.newFixedThreadPool(threadNum);
        ArrayList<Future<ILPsearch.Entry>> taskList = new ArrayList<>(datasetReader.getUsefulSpectraNum() + 10);
        Connection sqlConSpecCoder = DriverManager.getConnection(sqlPath);
        Statement sqlState = sqlConSpecCoder.createStatement();
        ResultSet sqlResSet = sqlState.executeQuery("SELECT scanName, scanNum, precursorCharge, precursorMass, mgfTitle FROM spectraTable");

        ReentrantLock lock = new ReentrantLock();
        int submitNum = 0;

        while (sqlResSet.next()) {
            String scanName = sqlResSet.getString("scanName");
            int scanNum = sqlResSet.getInt("scanNum");
            int precursorCharge = sqlResSet.getInt("precursorCharge");
            double precursorMass = sqlResSet.getDouble("precursorMass");
            if (java.lang.management.ManagementFactory.getRuntimeMXBean().getInputArguments().toString().indexOf("jdwp") >= 0){
                if (!debugScanNum.contains(scanNum)) {
                    continue;
                }
            }

            int fileId = fileNameIdMap.get( scanName.split("\\.")[0] );
            submitNum++;
            taskList.add(threadPool.submit(new ILPsearch(scanNum, buildIndex, massTool, ms2Tolerance, ms1Tolerance, inferPTM.getMinPtmMass(), inferPTM.getMaxPtmMass()
                    , spectraParserArray[fileId], minClear, maxClear, lock, scanName, precursorCharge, precursorMass, specProcessor, scanNum_truth_map.get(scanNum) )));
        }
        System.out.println("totalSubmit in SpecCoder, "+ submitNum);
        sqlResSet.close();
        sqlState.close();

        Map<Integer, Peptide> scanNumPeptideMap = new HashMap<>();
        Map<Integer, String> scanNumPepCandiStrMap = new HashMap<>();
        Map<VarPtm, Integer> varPtmCountMap = new HashMap<>();
        int lastProgress = 0;
        int totalCount = taskList.size();
        int count = 0;
        List<String> recStringList = new ArrayList<>(submitNum);
        while (count < totalCount) {
            // record search results and delete finished ones.
            List<Future<ILPsearch.Entry>> toBeDeleteTaskList = new ArrayList<>(totalCount - count);
            for (Future<ILPsearch.Entry> task : taskList) {
                if (task.isDone()) {
                    if (task.get() != null ) {
                        ILPsearch.Entry entry = task.get();
                        for (VarPtm varPtm : entry.varPtmList){
                            if (varPtmCountMap.containsKey(varPtm)) {
                                varPtmCountMap.put(varPtm, varPtmCountMap.get(varPtm)+1);
                            }else{
                                varPtmCountMap.put(varPtm, 1);
                            }
                        }
                        scanNumPeptideMap.put(entry.scanNum, entry.topPeptide);
                        scanNumPepCandiStrMap.put(entry.scanNum, entry.peptideSet);
                        recStringList.add(entry.recString);
                    }
                    toBeDeleteTaskList.add(task);
                } else if (task.isCancelled()) {
                    toBeDeleteTaskList.add(task);
                }
            }
            count += toBeDeleteTaskList.size();
            taskList.removeAll(toBeDeleteTaskList);
            taskList.trimToSize();

            int progress = count * 20 / totalCount;
            if (progress != lastProgress) {
                logger.info("Solving ILPs {}%...", progress * 5);
                lastProgress = progress;
            }

            if (count == totalCount) {
                break;
            }
            Thread.sleep(6000);
        }
        // shutdown threads.
        threadPool.shutdown();
        if (!threadPool.awaitTermination(60, TimeUnit.SECONDS)) {
            threadPool.shutdownNow();
            if (!threadPool.awaitTermination(60, TimeUnit.SECONDS))
                throw new Exception("Pool did not terminate");
        }
        if (lock.isLocked()) {
            lock.unlock();
        }

        FileWriter fileWriter;
        fileWriter = new FileWriter(outputDir+"ILP238recString"+".csv");
        BufferedWriter writer = new BufferedWriter(fileWriter);

        writer.write("scanNum,tagLen,truthPepLen,numVar\n");
        for (String recStr : recStringList) {
            writer.write(recStr);
        }
        writer.close();





        postProcessing(outputDir, massTool, varPtmCountMap, scanNumPeptideMap, scanNumPepCandiStrMap);
        logger.info("Saving results...");
    }

    class ScanRes {
        public double expMass;
        public int scanNum;
        public List<CandiScore> peptideInfoScoreList;
        public String scanName;
        ScanRes( int scanNum, List<CandiScore> peptideInfoScoreList, double expMass){
            this.scanNum = scanNum;
            this.peptideInfoScoreList = peptideInfoScoreList;
            this.expMass = expMass;
        }
    }
    class CandiScore implements Comparable<CandiScore>{
        public String ptmContainingSeq;
        public double pepScore;
        public double protScore = 0;
        public double varPtmTotalScore = 0;
        CandiScore(double pepScore, String ptmContainingSeq) {
            this.pepScore = pepScore;
            this.ptmContainingSeq = ptmContainingSeq;
        }

        public void setVarPtmTotalScore(Map<String, Double> varPtmStrScoreRefMap) { //varPtmScoreRefMap is like K(14.016)->1.2
            int startI = -1;
            double varPtmTotalScore = 0;
            int n_PtmOnPep = 0;
            for (int anyI = 0; anyI < ptmContainingSeq.length(); anyI++) {
                char thisChar = ptmContainingSeq.charAt(anyI);
                if (thisChar == '(') {
                    n_PtmOnPep++;
                    startI = anyI-1;
                } else if (thisChar == ')') {
                    String thisVarPtmStr = ptmContainingSeq.substring(startI, anyI + 1);
                    if (varPtmStrScoreRefMap.containsKey(thisVarPtmStr)) {
                        varPtmTotalScore += varPtmStrScoreRefMap.get(thisVarPtmStr);
                    }
                }
            }
            this.varPtmTotalScore = varPtmTotalScore/n_PtmOnPep;
        }
        public int compareTo(CandiScore o2) {
            if (this.protScore < o2.protScore) {
                return -1;
            } else if (this.protScore > o2.protScore) {
                return 1;
            } else {
                if (this.varPtmTotalScore < o2.varPtmTotalScore) {
                    return -1;
                } else if (this.varPtmTotalScore > o2.varPtmTotalScore) {
                    return 1;
                } else {
                    if (this.pepScore < o2.pepScore) {
                        return -1;
                    } else if (this.pepScore > o2.pepScore) {
                        return 1;
                    }
                }
            }
            return 0;
        }
    }
    private void postProcessing(String outputDir,  MassTool massTool, Map<VarPtm, Integer> varPtmCountMap, Map<Integer, Peptide> scanNumPeptideMap,Map<Integer, String> scanNumPepCandiStrMap) throws IOException {
        //preprocess varPtmCountMap
        Map<String, Double> varPtmRefScoreMap = new HashMap<>();
        for (VarPtm varPtm : varPtmCountMap.keySet()) {
            if (varPtmCountMap.get(varPtm) == 1) {
                continue;
            }
            varPtmRefScoreMap.put(varPtm.site+"("+InferPTM.df3.format(varPtm.mass)+")", Math.sqrt(varPtmCountMap.get(varPtm)));
        }
        List<Map.Entry<String, Double>> testList = new ArrayList<>(varPtmRefScoreMap.entrySet());
        Collections.sort(testList, Map.Entry.comparingByValue(Comparator.reverseOrder()));
        int j = 0;
        for (Map.Entry<String, Double> entry : testList) {
            String varPtmStr = entry.getKey();
            if (j>18) {
                varPtmRefScoreMap.remove(varPtmStr);
            }
            j++;
        }

        for (String varPtmStr : varPtmRefScoreMap.keySet()) {
            System.out.println("after,"+varPtmStr + "," + InferPTM.df3.format(varPtmRefScoreMap.get(varPtmStr)));
        }
        //collect data

        List<ScanRes> scanResList = new ArrayList<>();
        for (int scanNum : scanNumPeptideMap.keySet()) {
            Peptide topPeptide = scanNumPeptideMap.get(scanNum);
            String peptideSet = scanNumPepCandiStrMap.get(scanNum);
            String[] candiSetStr = peptideSet.split(",");
            int numPep = candiSetStr.length/3;

            List<CandiScore> candiScoreList = new ArrayList<>();
            for (int i = 0; i < numPep; i++) {
                String ptmContainingSeq = candiSetStr[3*i+0];
                double thisScore = Double.valueOf(candiSetStr[3*i+1]);
                CandiScore candiScore = new CandiScore(thisScore, ptmContainingSeq);
                candiScore.setVarPtmTotalScore(varPtmRefScoreMap);
                candiScoreList.add(candiScore); //peptideInfo and their score
            }
            Collections.sort(candiScoreList, Comparator.comparing(o -> o.pepScore, Comparator.reverseOrder()));
            double topPepScore = candiScoreList.get(0).pepScore;
            Iterator<CandiScore> iter = candiScoreList.iterator();
            while (iter.hasNext()) {
                CandiScore candiScore = iter.next();
                if (candiScore.pepScore < 0.85 * topPepScore) {
                    iter.remove();
                }
            }
            scanResList.add(new ScanRes(scanNum, candiScoreList, topPeptide.precursorMass));
        }

        for (ScanRes scanRes : scanResList) {
            Collections.sort(scanRes.peptideInfoScoreList, Comparator.reverseOrder()); // rank candidates using peptide score
        }
        Collections.sort(scanResList, Comparator.comparing(o -> o.peptideInfoScoreList.get(0).pepScore, Comparator.reverseOrder()));


        DecimalFormat df= new  DecimalFormat( ".00000" );
        List<Pair<Double, String>> finalExcelList = new ArrayList<>(scanResList.size());
        for (ScanRes scanRes : scanResList) {
            if (debugScanNum.contains(scanRes.scanNum)) {
                int a = 1;
            }
            List<CandiScore> candiScoreList = scanRes.peptideInfoScoreList;
            Peptide topPep = scanNumPeptideMap.get(scanRes.scanNum);
            CandiScore topCandi = candiScoreList.get(0);
            double topPepScore = topCandi.pepScore;

            double theoMass = massTool.calResidueMass(topCandi.ptmContainingSeq) + massTool.H2O;
            double massDiff = getMassDiff(scanRes.expMass, theoMass, MassTool.C13_DIFF);
            double ppm = Math.abs(massDiff * 1e6 / theoMass);
            StringBuilder fullPtmSeq = new StringBuilder();
//            StringBuilder fullFreeSeq = new StringBuilder();
            for (CandiScore candi : candiScoreList) {
                if (Math.abs(candi.pepScore-topPepScore) < 0.001) {
                    fullPtmSeq.append(candi.ptmContainingSeq).append(";");
//                    fullFreeSeq.append(topPep.getFreeSeq()).append(";");
                }
            }
            String finalStr = String.format(Locale.US, "%d,%s,%s,%f,%f,%f\n"
                    , scanRes.scanNum,  fullPtmSeq,  df.format(topCandi.pepScore)
                    ,  ppm, theoMass, scanRes.expMass
            );

            finalExcelList.add(new Pair(topCandi.pepScore, finalStr));
        }

        Collections.sort(finalExcelList, Comparator.comparing(o -> o.getFirst(), Comparator.reverseOrder()));
        FileWriter fileWriter;
        try {
            fileWriter = new FileWriter(outputDir+"ILP238"+".csv");
        } catch (Exception e){
            fileWriter = new FileWriter(outputDir+"ILP238"+ df.format(Math.random())+".csv");
        }
        BufferedWriter writer = new BufferedWriter(fileWriter);

        writer.write("scanNum,peptide,pepScore,ppm,theoMass,expMass\n");
        for (Pair<Double, String> pair : finalExcelList) {
            writer.write(pair.getSecond());
        }
        writer.close();
    }

    public static double getMassDiff(double expMass, double theoMass, double C13Diff) {
        double massDiff1 = expMass - theoMass;
        double massDiff2 = expMass - theoMass - C13Diff;
        double massDiff3 = expMass - theoMass - 2 * C13Diff;
        double absMassDiff1 = Math.abs(massDiff1);
        double absMassDiff2 = Math.abs(massDiff2);
        double absMassDiff3 = Math.abs(massDiff3);

        if ((absMassDiff1 <= absMassDiff2) && (absMassDiff1 <= absMassDiff2)) {
            return massDiff1;
        } else if ((absMassDiff2 <= absMassDiff1) && (absMassDiff2 <= absMassDiff3)) {
            return massDiff2;
        } else {
            return massDiff3;
        }
    }

}
