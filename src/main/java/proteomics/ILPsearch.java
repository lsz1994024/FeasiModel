package proteomics;

import ProteomicsLibrary.MassTool;
import ProteomicsLibrary.Score;
import ProteomicsLibrary.SpecProcessor;
import ProteomicsLibrary.Types.SparseBooleanVector;
import ProteomicsLibrary.Types.SparseVector;
import com.google.common.collect.Sets;
import gurobi.GRB;
import gurobi.GRBEnv;
import org.apache.commons.math3.util.Pair;
import proteomics.FM.FMRes;
import proteomics.Index.BuildIndex;
import proteomics.PTM.InferPTM;
import proteomics.Segment.InferSegment;
import proteomics.Spectrum.DatasetReader;
import proteomics.Types.*;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;

import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.locks.ReentrantLock;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static proteomics.PIPI.*;
import static proteomics.PTM.InferPTM.*;
import static proteomics.Segment.InferSegment.*;

public final class ILPsearch implements Callable<ILPsearch.Entry> {
    private static final int candisNum = 20;
    private final BuildIndex buildIndex;
    private final MassTool massTool;
    private final double ms1Tolerance;
    private final double minPtmMass;
    private final double maxPtmMass;
    private final JMzReader spectraParser;
    private final double minClear;
    private final double maxClear;
    private final ReentrantLock lock;
    private final String scanName;
    private final int precursorCharge;
    private final double precursorMass;
    private final SpecProcessor specProcessor;
    private final int scanNum;
    private final double ms2Tolerance;
    private final InferPTM inferPTM;
    private List<String> truthStrList = null;
    public ILPsearch(int scanNum, BuildIndex buildIndex, MassTool massTool, double ms2Tolerance, double ms1Tolerance, double minPtmMass, double maxPtmMass
            , JMzReader spectraParser, double minClear, double maxClear, ReentrantLock lock, String scanName, int precursorCharge, double precursorMass
            , SpecProcessor specProcessor, List<String> truthStrList) {

        this.buildIndex = buildIndex;
        this.massTool = massTool;
        this.ms1Tolerance = ms1Tolerance;
        this.minPtmMass = minPtmMass;
        this.maxPtmMass = maxPtmMass;
        this.truthStrList = truthStrList;
        this.spectraParser = spectraParser;
        this.minClear = minClear;
        this.maxClear = maxClear;
        this.lock = lock;
        this.scanName = scanName;
        this.precursorCharge = precursorCharge;
        this.precursorMass = precursorMass;
        this.specProcessor = specProcessor;
        this.scanNum = scanNum;
        this.ms2Tolerance = ms2Tolerance;
        this.inferPTM = buildIndex.getInferPTM();
    }

    @Override
    public Entry call() throws Exception {
        Map<Double, Double> rawPLMap;
        try {
            lock.lock();
            rawPLMap = spectraParser.getSpectrumById(scanName.split("\\.")[1]).getPeakList();
        } finally {
            lock.unlock();
        }

        double ms1TolAbs = Double.parseDouble(InferPTM.df3.format(precursorMass*ms1Tolerance/1000000));

        TreeMap<Double, Double> plMap = specProcessor.preSpectrumTopNStyleWithChargeLimit(rawPLMap, precursorMass, precursorCharge, minClear, maxClear, DatasetReader.topN, ms2Tolerance);

        if (plMap.isEmpty()) return null;
        // Coding
        SparseVector expProcessedPL = specProcessor.digitizePL(plMap);
        InferSegment inferSegment = buildIndex.getInferSegment();
        TreeMap<Double, Double> finalPlMap = inferSegment.addVirtualPeaks(precursorMass, plMap);

//        try{
//            inferSegment.getOneCorrectTag(truthStrList);
//        } catch (Exception e) {
//            System.out.println(scanNum+", wrong");
//        }
        ExpTag refTag = inferSegment.getOneCorrectTag(truthStrList);
        if (debugScanNum.contains(this.scanNum)) {
            int a = 1;
        }
        if (refTag == null)  return null;
        Recorder rec = new Recorder(scanNum, refTag.size());

        double totalMass = precursorMass + 2 * MassTool.PROTON;


        Map<String, PeptideInfo> peptideInfoMap = new HashMap<>(50000);

        GRBEnv env = new GRBEnv(true);
        env.set(GRB.IntParam.OutputFlag,0);
        env.set(GRB.IntParam.LogToConsole, 0);
        env.start();

        Map<String, String> protSeqMap = buildIndex.protSeqMap;
        TreeSet<Peptide> peptideTreeSet = new TreeSet<>(Collections.reverseOrder());
        List<OccGroup> occGList = new ArrayList<>(100);
        getOccGroupList(scanNum , refTag, occGList);
        occGList.sort(Comparator.comparing(o->o.totalScore, Comparator.reverseOrder()));

        String truthFreeSeq = truthStrList.get(0).replace('I', 'L');

        int count = 0;
        boolean truthIncluded = false;
        outLoop:
        for (OccGroup occG : occGList) {
            Pair<ExpTag, Integer> tagRelPos = occG.tagRelPos;
            ExpTag tag = tagRelPos.getFirst();
            int relPos = tagRelPos.getSecond();
            String protId = occG.protId;
            if (relPos+truthFreeSeq.length() < protSeqMap.get(protId).length()
                    && ( protId.contains("Packet_Kmod:")
                    && protSeqMap.get(protId).substring(relPos, relPos+truthFreeSeq.length()).contentEquals(truthFreeSeq) )
            ) {
                truthIncluded = true;
            }
            count++;
            if (count > 20) {
                if (! truthIncluded) {
                    System.out.println(scanNum +", truth prot and pos not included.");
                }
                break outLoop;
            }
            addCandisWithMultiMaxPeakMILP(scanNum, protId, tag, relPos,  ms1TolAbs, peptideTreeSet, peptideInfoMap, expProcessedPL, plMap, env, rec);
        }
        if (debugScanNum.contains(scanNum)){
            int a = 1;
        }
        if (peptideTreeSet.isEmpty()) {
            env.dispose();
            return null;
        }

        Iterator<Peptide> iterator = peptideTreeSet.iterator();
        while (iterator.hasNext()) {
            Peptide peptide = iterator.next();
            int numMC = getNumOfMCFromStr(peptide.getFreeSeq().substring(0, peptide.getFreeSeq().length()-1));
            if (numMC > 2) {
                iterator.remove();
            }
        }

        List<Peptide> pepList = new ArrayList<>(peptideTreeSet);
        Set<Integer> pepIdsToRemove = new HashSet<>();
        for (int j = pepList.size()-1; j >= 1; j--) {
            if (pepIdsToRemove.contains(j)) continue;
            for (int i = 0; i < j; i++) {
                if (pepIdsToRemove.contains(i)) continue;
                if (isHomo(pepList.get(i), pepList.get(j), peptideInfoMap) || pepList.get(j).getScore() > 0.6*pepList.get(i).getScore()) {
                    int iPriority = pepList.get(i).getPriority();
                    int jPriority = pepList.get(j).getPriority();
                    if (iPriority < jPriority) {
                        pepIdsToRemove.add(i);
                    } else if (iPriority > jPriority) {
                        pepIdsToRemove.add(j);
                    } else {// iPriority == jPriority

                        if (onlyDifferUnsettledPtm(pepList.get(i), pepList.get(j))) {
                            pepIdsToRemove.add(j);
                        }
                    }
                }
            }
        }
        List<Peptide> newPepList = new ArrayList<>();
        for (int id = 0; id < pepList.size(); id++){
            if (pepIdsToRemove.contains(id)) continue;
            newPepList.add(pepList.get(id));
        }
        if (newPepList.isEmpty()) {
            env.dispose();
            return null;
        }

        Peptide[] peptideArray = newPepList.toArray(new Peptide[0]);
        Peptide topPep = peptideArray[0];

        Entry entry;
        String pepSetString = "";
        for (Peptide peptide : peptideArray){
            PeptideInfo peptideInfo = peptideInfoMap.get(peptide.getFreeSeq());
            pepSetString += peptide.getVarPtmContainingSeqNow() + "," + peptide.getScore() + "," + String.join("_", peptideInfo.protIdSet) +",";
        }

        double deltaLCn = 1; // L means the last?
        if (peptideArray.length > candisNum - 1) {
            deltaLCn = (peptideArray[0].getScore() - peptideArray[candisNum - 1].getScore()) / peptideArray[0].getScore();
        }
        double deltaCn = 1;
        if (peptideArray.length > 1) {
            for(int i = 0; i < peptideArray.length; i++) {
                if (peptideArray[i].getScore() != peptideArray[0].getScore()){
                    deltaCn = (peptideArray[0].getScore() - peptideArray[i].getScore()) / peptideArray[0].getScore();
                    break;
                }
            }
        }
        String otherPtmPatterns = "-";
        entry = new ILPsearch.Entry(
                scanNum, scanName
                , precursorCharge, precursorMass, buildIndex.getLabelling(), topPep.getPtmContainingSeq(buildIndex.returnFixModMap())
                , topPep.getTheoMass(), topPep.isDecoy() ? 1 : 0, topPep.getScore(), deltaLCn, deltaCn
                , topPep.getMatchedPeakNum(), topPep.getIonFrac(), topPep.getMatchedHighestIntensityFrac()
                , topPep.getExplainedAaFrac(), otherPtmPatterns, topPep.getaScore(), ""
                , pepSetString.substring(0, pepSetString.length()-1)
        );
        entry.topPeptide = topPep;
        entry.topPeptide.precursorMass = precursorMass;
        entry.varPtmList.addAll(topPep.posVarPtmResMap.values());
        entry.recString = scanNum + "," + rec.tagLen + ","+ truthStrList.get(0).length() +","+ rec.x_num + "\n";
        for (Peptide peptide : peptideArray) {
            entry.peptideInfoMapForRef.put(peptide.getFreeSeq(), peptideInfoMap.get(peptide.getFreeSeq()));
        }
        env.dispose();
        return entry;
    }

    final class Recorder {
        public int x_num = 0;
        public int scanNum = 0;
        public int tagLen = 0;
        public Recorder(int scanNum, int tagLen) {
            this.scanNum = scanNum;
            this.tagLen = tagLen;
        }
    }
    final class OccGroup {
        public Pair<ExpTag, Integer> tagRelPos;
        final public String protId;
        public double totalScore;
        public OccGroup(ExpTag initialTag, String protId, int tagPosInProt) {
            this.tagRelPos = new Pair<>(initialTag, tagPosInProt);
            this.protId = protId;
            this.totalScore = initialTag.getTotalIntensity();
        }
    }
    private void getProtCandidates(int scanNum, ExpTag refTag, List<OccGroup> occGroupList, Map<String, String> protSeqMap) {
        String tagStr = refTag.getFreeAaString();
        for (String protId : protSeqMap.keySet()) {
            String protSeq = protSeqMap.get(protId);
            int tagPosInProt = protSeq.indexOf(tagStr);
            while (tagPosInProt >= 0) {
                OccGroup occG = new OccGroup(refTag, protId, tagPosInProt);
                occGroupList.add(occG);
                tagPosInProt = protSeq.indexOf(tagStr, tagPosInProt+1);
            }
        }
    }
    private boolean isHomo(Peptide p1, Peptide p2, Map<String, PeptideInfo> peptideInfoMap) {
//        int a = 1;
        Set<String> temp1 = peptideInfoMap.get(p1.getFreeSeq()).protIdSet;
        Set<String> temp2 = peptideInfoMap.get(p2.getFreeSeq()).protIdSet;

        Set<String> set = new HashSet<>(temp1);
        set.retainAll(temp2);
        SparseBooleanVector sbv1 = buildIndex.inferSegment.generateSegmentBooleanVector(p1.getFreeSeq());
        SparseBooleanVector sbv2 = buildIndex.inferSegment.generateSegmentBooleanVector(p2.getFreeSeq());
        return sbv1.dot(sbv2) > 0.3*Math.min(p1.getFreeSeq().length(), p2.getFreeSeq().length());
    }

    private int getOccGroupList(int scanNum, ExpTag tag, List<OccGroup> occGroupList) {
        char[] tagChar;
        int numRes;
        int solCount = 0;
        FMRes fmRes;
        String protId;
        tagChar = tag.getFreeAaString().toCharArray();
        fmRes = buildIndex.fmIndexFull.fmSearch(tagChar);
        solCount += fmRes.ep-fmRes.sp+1;

        int absTagPos;
        int dotIndex;
        int lPos;
        for (int ii = fmRes.sp; ii <= fmRes.ep; ii++) {
            absTagPos = buildIndex.fmIndexFull.SA[ii];
            dotIndex = -2-Arrays.binarySearch(buildIndex.dotPosArrFull, absTagPos);
            protId = buildIndex.posProtMapFull.get(dotIndex);
            lPos = absTagPos - buildIndex.dotPosArrFull[dotIndex] - 1;

            OccGroup occG = new OccGroup(tag, protId, lPos);
            occGroupList.add(occG);
        }
        return solCount;
    }

    private boolean onlyDifferUnsettledPtm(Peptide p1, Peptide p2) {
        SparseBooleanVector sbv1 = buildIndex.inferSegment.generateSegmentBooleanVector(p1.getFreeSeq());
        SparseBooleanVector sbv2 = buildIndex.inferSegment.generateSegmentBooleanVector(p2.getFreeSeq());
        if (sbv1.dot(sbv2) < 0.3*Math.min(p1.getFreeSeq().length(), p2.getFreeSeq().length())){
            return false;
        }

        if (p1.posVarPtmResMap.size() != p2.posVarPtmResMap.size()) {
            return false;
        }
        if (!p1.posVarPtmResMap.keySet().containsAll(p2.posVarPtmResMap.keySet())){
            return false;
        }
        byte n_SameMass = 0;
        for (int pos : p2.posVarPtmResMap.keySet()) {
            if (p1.posVarPtmResMap.get(pos).mass == p2.posVarPtmResMap.get(pos).mass) {
                n_SameMass++;
            } else if (Math.abs(p1.posVarPtmResMap.get(pos).mass - p2.posVarPtmResMap.get(pos).mass) < 0.02) {
                if (p1.posVarPtmResMap.get(pos).classification.contains("PIPI_unsettled") || p2.posVarPtmResMap.get(pos).classification.contains("PIPI_unsettled")) {
                    n_SameMass++;
                }
            }
        }
        return n_SameMass == p1.posVarPtmResMap.size();
    }

    private int getNumOfMCFromStr(String tagStr) {
        String str2 = tagStr.replaceAll("[KR]","");
        int numMC = tagStr.length()-str2.length();
        return numMC;
    }
    private double addCandisWithMultiMaxPeakMILP(int scanNum, String protId, ExpTag refTag, int tagPosInProt, double ms1TolAbs, TreeSet<Peptide> resPepTreeSet
            , Map<String, PeptideInfo> peptideInfoMap, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, GRBEnv env, Recorder rec) {

        int tagNum;
        String protSeq = buildIndex.protSeqMap.get(protId);
        int protLen = protSeq.length();
        if (tagPosInProt < 0 || (tagPosInProt == 0 && refTag.isNorC != N_TAG)) {
            return 0;
        }
        if (tagPosInProt >= protLen
                || tagPosInProt + refTag.size() == protLen && refTag.isNorC != C_TAG) {
            return 0;
        }

        double maxScore = 0;
        int remainMC = massTool.missedCleavage - getNumOfMCFromStr(refTag.getFreeAaString());
        if (refTag.isNorC == C_TAG && isKR(refTag.getFreeAaString().charAt(refTag.size()-1))) {
            remainMC++; // the last KR at pepC does not count as MC
        }

        List<TreeSet<Peptide>> segResList = new ArrayList<>();
        boolean solved;

        //N
        if (Math.abs(refTag.getHeadLocation()-MassTool.PROTON) < ms2Tolerance) {//rTag.isNorC == N_TAG
            solved = true;
        } else {
            TreeSet<Peptide> nModPepsSet = new TreeSet<>(Comparator.reverseOrder());
            solved = solveGapN(scanNum, refTag, tagPosInProt, protSeq, remainMC, ms1TolAbs, expProcessedPL, plMap, env, nModPepsSet, rec);
            if (solved) segResList.add(nModPepsSet);
        }
        if (! solved) return 0;

        //Tag
        if (solved) {
            segResList.add(getPepFromTag(refTag));
        } else {
            return 0;
        }

        //C
        if (Math.abs(refTag.getTailLocation()-precursorMass+ massTool.H2O-MassTool.PROTON) < ms1TolAbs+ms2Tolerance) {// lTag.isNorC == C_TAG
            solved = true;
        } else {
            TreeSet<Peptide> cModPepsSet = new TreeSet<>(Comparator.reverseOrder());
            solved = solveGapC(scanNum, refTag, tagPosInProt, protId, protSeq, remainMC,
                    ms1TolAbs, expProcessedPL, plMap, env, cModPepsSet, rec);
            if (solved) segResList.add(cModPepsSet);
        }
        if (!solved) { // todo makeup solution
            return 0;
        }
        if ( ! segResList.isEmpty()) {
            collectResult(segResList, resPepTreeSet, expProcessedPL, plMap,  peptideInfoMap, protId, protSeq);
        }
        return maxScore;
    }

    private void collectResult(List<TreeSet<Peptide>> segResList, TreeSet<Peptide> resPepTreeSet, SparseVector expProcessedPL, TreeMap<Double,Double> plMap,
                               Map<String, PeptideInfo> peptideInfoMap, String protId, String protSeq) {
        List<Set<Integer>> idSetList = new ArrayList<>(segResList.size());
        for (int segNum = 0; segNum < segResList.size(); segNum++) {
            Set<Integer> idSet = new HashSet<>();
            for (int id = 0; id < segResList.get(segNum).size(); ++id) {
                idSet.add(id);
            }
            idSetList.add(idSet);
        }

        for (List<Integer> resIdList : Sets.cartesianProduct(idSetList)) {
            StringBuilder pepSeqSB = new StringBuilder();
            PosMassMap posMassMap = new PosMassMap();
            TreeMap<Integer, VarPtm> posVarPtmResMap = new TreeMap<>();

            int curRelPos = 0;
            for (int i = 0; i < resIdList.size(); ++i) {
                int id = resIdList.get(i);
                Peptide[] segArr = new Peptide[segResList.get(i).size()];
                segResList.get(i).toArray(segArr);
                Peptide curPepSeg = segArr[id];
                pepSeqSB.append(curPepSeg.getFreeSeq());

                if (! curPepSeg.posVarPtmResMap.isEmpty()) {
                    for (int j : curPepSeg.posVarPtmResMap.keySet()) {
                        posMassMap.put(j+curRelPos, curPepSeg.posVarPtmResMap.get(j).mass);
                        posVarPtmResMap.put(j+curRelPos, curPepSeg.posVarPtmResMap.get(j));
                    }
                }
                curRelPos += segArr[id].getFreeSeq().length();
            }
            Peptide resPep = new Peptide(pepSeqSB.toString(), false, massTool);
            if (! posMassMap.isEmpty()) {
                resPep.setVarPTM(posMassMap);
                resPep.posVarPtmResMap.putAll(posVarPtmResMap);
            }
            double calScore = massTool.buildVectorAndCalXCorr(resPep.getIonMatrixNow(), 1, expProcessedPL, resPep.matchedBions, resPep.matchedYions);//todo decide the penalty

            resPep.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, resPep.getIonMatrixNow(), ms2Tolerance));
            resPep.setScore(calScore*(1 - resPep.posVarPtmResMap.size() * 0.02));

            updatePeptideTreeSet(resPep, resPepTreeSet, peptideInfoMap, protId, protSeq, protSeq.indexOf(pepSeqSB.toString()), protSeq.indexOf(pepSeqSB.toString())+pepSeqSB.length()-1);
        }
    }

    private TreeSet<Peptide> getPepFromTag(ExpTag tag) {
        String freeSeq = tag.getFreeAaString();
        String modSeq  = tag.getPtmAaString();
        Peptide peptide = new Peptide(freeSeq, false, massTool);
        TreeSet<Peptide> modPepPool = new TreeSet<>();
        modPepPool.add(peptide);
        if (! freeSeq.contentEquals(modSeq)) {
            PosMassMap posMassMap = new PosMassMap();
            int idOfAa = -1;
            for (char aaChar : tag.getPtmAaString().toCharArray()) {
                if (Character.isUpperCase(aaChar)) {
                    idOfAa += 1;
                } else {
                    posMassMap.put(idOfAa, massTool.labelVarPtmMap.get(aaChar).mass);
                    peptide.posVarPtmResMap.put(idOfAa, massTool.labelVarPtmMap.get(aaChar));
                }
            }
        }
        return modPepPool;
    }
    private boolean solveGapC(int scanNum, ExpTag tag, int tagPosInProt, String protId, String protSeq, int remainMC, double ms1TolAbs,
                              SparseVector unUsedExpProcessedPL, TreeMap<Double, Double> unUsedPlMap, GRBEnv env, TreeSet<Peptide> cModPepsSet, Recorder rec) {

        double tagCMass = tag.getTailLocation() + massTool.H2O - MassTool.PROTON; // +massTool.H2O-MassTool.PROTON  is to mimic the mass of the real neutral precursor mass
        if (tag.isNorC == C_TAG || Math.abs(tagCMass - precursorMass) < ms1TolAbs+ms2Tolerance) return true;// C tag, tagCMass must be exactly settle. Otherwise no valid cPos will be, then no valid n pos.
        // only when oriTag is not C oriTag can it extend to c
        double cCutMass = tag.getTailLocation() - MassTool.PROTON;
        int tagLen = tag.size();
        int protLen = protSeq.length();
        List<Integer> cPoses = new ArrayList<>();
        List<Integer> krPoses = new ArrayList<>();
        int mcInC = 0;
        boolean startRecord = false;
        char aaChar;
        boolean isKR;

        for (int cPos = tagPosInProt+tagLen; cPos < protLen; cPos++) {  //
            aaChar = protSeq.charAt(cPos);
            if (isX(aaChar))  break;
            isKR = isKR(aaChar);
            tagCMass += massTool.getMassTable().get(aaChar); //even when tag is fuzzy, tagcmass wont be disturbed
            if (mcInC > remainMC || tagCMass > precursorMass+maxPtmMass) break;
            if (tagCMass >= precursorMass-maxPtmMass) {
                if (startRecord) {
                    cPoses.add(cPos);
                } else {
                    startRecord = true;
                }
                if (isKR)  krPoses.add(cPos);

                if (Math.abs(tagCMass - precursorMass) <= ms1TolAbs+ms2Tolerance) {
                    String cPartSeq = protSeq.substring(tagPosInProt+tagLen, cPos+1);
                    storeCleanPartPeptides(cPartSeq, cCutMass, C_PART, unUsedExpProcessedPL, unUsedPlMap, cModPepsSet);
                    return true;
                }
            }
            if (isKR) mcInC++;
        }// finish collect all possible cposes and krPoses
        if ( !startRecord) return false;// cPoses is empty

        Map<Integer, Integer> posYIdMap = new HashMap<>();
        Map<Integer, Integer> yIdMaxAbsPosMap = new HashMap<>();
        int optStartPos = 1;
        int optEndPosP1 = 0;
        if ( ! krPoses.isEmpty()) {
            optEndPosP1 = krPoses.get(krPoses.size()-1) + 1;  //max kr
            optStartPos = krPoses.get(0) + 1;  // good trick//min kr
            if ( ! cPoses.isEmpty() && cPoses.get(cPoses.size()-1) == protLen - 1) {
                optEndPosP1 = protLen;
            }
            int yId = 0;
            for (int pos = optStartPos; pos < optEndPosP1; pos++) {
                posYIdMap.put(pos, yId); // using rel pos in part seq i.e. start from 0 not the prot pos
                int oldMaxPos = yIdMaxAbsPosMap.getOrDefault(yId, -99);
                yIdMaxAbsPosMap.put(yId, pos > oldMaxPos ? pos : oldMaxPos);
                if (krPoses.contains(pos)) {
                    yId++;
                }
            }
        } else {
            if (cPoses.isEmpty()) {
                return false;
            } else if (cPoses.get(cPoses.size()-1) == protLen-1) {
                optEndPosP1 = protLen;
                optStartPos = protLen;  // good trick
            }
        }

        if (optStartPos <= optEndPosP1) {
            int cPartStartPos = tagPosInProt + tagLen;
            String cPartSeq = protSeq.substring(cPartStartPos, optEndPosP1);
            int cPartSeqLen = cPartSeq.length();
            double flexiableMass = precursorMass - (tag.getTailLocation() + massTool.H2O-MassTool.PROTON) - massTool.calResidueMass(protSeq.substring(cPartStartPos, optStartPos));
            Map<Integer, TreeMap<Double, VarPtm>> absPos_MassVarPtm_Map = new HashMap<>(cPartSeqLen, 1);
            Map<Integer, List<Byte>> absPos_ptmPositions_Map = new HashMap<>(cPartSeqLen, 1);
            inferPTM.prepareInfoCTerm(scanNum, cPartSeq, absPos_MassVarPtm_Map,  yIdMaxAbsPosMap, optEndPosP1 == protLen, optStartPos, cPartStartPos, absPos_ptmPositions_Map);

            if (cPartSeqLen == 1) {
                findPtmOnOneAa(cModPepsSet, flexiableMass, absPos_MassVarPtm_Map, cPartSeq, cPartStartPos, ms1TolAbs);
            } else {
                Map<Integer, Set<Double>> oneTimeMassGroups = new HashMap<>(cPartSeqLen);
                Map<Set<Integer>, Set<Double>> posComb_multiMassSet_Map = new HashMap<>(200);
                Map<Integer, Map<Double, VarPtm>> pos_MassVarPtm_Map = new HashMap<>(cPartSeqLen, 1);
                Map<Double, Set<Integer>> allMassAllPosesMap = new HashMap<>(cPartSeqLen * 30, 1); //Set<pos>


                inferPTM.prepareInfoCTerm(scanNum, cPartSeq, oneTimeMassGroups, posComb_multiMassSet_Map,
                        pos_MassVarPtm_Map, allMassAllPosesMap, yIdMaxAbsPosMap, optEndPosP1 == protLen, optStartPos, cPartStartPos);

                Set<Pair<Integer, Map<Double, Integer>>> resList = new HashSet<>(100);
                rec.x_num += inferPTM.findBestPtmMIPExtC(scanNum, env, allMassAllPosesMap, flexiableMass, (cPartStartPos),
                        cPartSeq, ms1TolAbs, oneTimeMassGroups, posComb_multiMassSet_Map, posYIdMap, pos_MassVarPtm_Map, resList);

                inferPTM.getFeasibleMassPosMapC(scanNum, resList, unUsedPlMap, cPartSeq, cCutMass, C_PART,
                        unUsedExpProcessedPL, false, allMassAllPosesMap, pos_MassVarPtm_Map, cModPepsSet, cPartStartPos, yIdMaxAbsPosMap, optStartPos);
            }
            if (yIdMaxAbsPosMap.isEmpty()) {
                inferPTM.findPosssible1Ptm(scanNum, cPartSeq, absPos_MassVarPtm_Map, cPartStartPos, cModPepsSet, flexiableMass,ms1TolAbs);
            }
        }

        if (cModPepsSet.isEmpty()) return false;

        return true;
    }


    private boolean solveGapN(int scanNum, ExpTag tag, int tagPosInProt, String protSeq, int remainMC, double ms1TolAbs,
                              SparseVector unUsedExpProcessedPL, TreeMap<Double, Double> unUsedPlMap, GRBEnv env, TreeSet<Peptide> nModPepsSet, Recorder rec) {
        double nDeltaMass = tag.getHeadLocation() - MassTool.PROTON; //n term delta mass should be independent of cDeltaMass
        double nCutMass = precursorMass + MassTool.PROTON - tag.getHeadLocation() - massTool.H2O;
        if (tag.isNorC == N_TAG || Math.abs(tag.getHeadLocation() - MassTool.PROTON) < ms1TolAbs)  return true;

        int mcInN = 0;
        char aaChar;

        //1 find fixedStartPos=optEndPosP1, optStartPos, only by mass
        int fixedStartPos = -1;
        int optStartPos = -1;
        for (int nPos = tagPosInProt-1; nPos >= 0; nPos--) {
            aaChar = protSeq.charAt(nPos);
            if (isX(aaChar)) break;
            nDeltaMass -= massTool.getMassTable().get(aaChar);
            if (Math.abs(nDeltaMass) <= ms1TolAbs) { // whenever this meets, store clean part seq
                String nPartSeq = protSeq.substring(nPos, tagPosInProt);
                storeCleanPartPeptides(nPartSeq, nCutMass, N_PART, unUsedExpProcessedPL, unUsedPlMap, nModPepsSet);
                return true;
            }
            if (nDeltaMass < maxPtmMass) { // only record the first time
                if (fixedStartPos == -1) {
                    fixedStartPos = nPos;
                    optStartPos = nPos;
                }
                if (nDeltaMass >= minPtmMass) {
                    optStartPos = nPos;
                } else {
                    break;
                }
            }
        }

        //2 update fixedStartPos, optStartPos and krPoses based on missedCleavageSite
        if (fixedStartPos == -1) return false; //not feasible
        for (int nPos = tagPosInProt-1; nPos >= fixedStartPos; nPos--) {
            aaChar = protSeq.charAt(nPos);
            if (isKR(aaChar)) {
                mcInN++;
                if (mcInN > remainMC) {
                    return false;
                }
            }
        }
        List<Integer> krPoses = new ArrayList<>();
        for (int nPos = fixedStartPos-1; nPos >= Math.max(0, optStartPos-1); nPos--) {
            aaChar = protSeq.charAt(nPos);
            if (isKR(aaChar)) {
                krPoses.add(nPos);
                mcInN++;
                if (mcInN > remainMC) {
                    optStartPos = nPos+1;
                    break;
                }
            }
        }

        //3 update fixedStartPos, optStartPos based on N-term specificity
        if (fixedStartPos >= 2) {
            if (krPoses.isEmpty()) {
                return false;
            } else {
                fixedStartPos = 1+krPoses.get(0); // the first krPos == the biggest krPose
                optStartPos = 1+krPoses.get(krPoses.size()-1); //the last krPos == the smallest krPose
            }
        }

        Map<Integer, Integer> posYIdMap = new HashMap<>();
        Map<Integer, Integer> yIdMinAbsPosMap = new HashMap<>();
        int yId = -1;

        if (fixedStartPos >= 2) {
            for (int pos = fixedStartPos - 1; pos >= optStartPos; pos--) {
                if (krPoses.contains(pos)) yId++;
                int tmpMinPos = yIdMinAbsPosMap.getOrDefault(yId, 999999);
                yIdMinAbsPosMap.put(yId, pos < tmpMinPos ? pos : tmpMinPos);
                posYIdMap.put(pos, yId);
            }
        } else {
            for (int pos = fixedStartPos - 1; pos >= optStartPos; pos--) {
                yId++;
                int tmpMinPos = yIdMinAbsPosMap.getOrDefault(yId, 999999);
                yIdMinAbsPosMap.put(yId, pos < tmpMinPos ? pos : tmpMinPos);
                posYIdMap.put(pos, yId);
            }
        }

        if (optStartPos <= fixedStartPos) {
            String nPartSeq = protSeq.substring(optStartPos, tagPosInProt);
            int nPartSeqLen = tagPosInProt - optStartPos;
            double flexiableMass = tag.getHeadLocation() - MassTool.PROTON - massTool.calResidueMass(protSeq.substring(fixedStartPos, tagPosInProt));

            Map<Integer, TreeMap<Double, VarPtm>> absPos_MassVarPtm_Map = new HashMap<>(nPartSeqLen, 1);
            Map<Integer, List<Byte>> absPos_ptmPositions_Map = new HashMap<>(nPartSeqLen, 1);

            inferPTM.prepareInfoNTerm(scanNum, nPartSeq,
                    absPos_MassVarPtm_Map, yIdMinAbsPosMap, optStartPos == 0, fixedStartPos, tagPosInProt, protSeq, absPos_ptmPositions_Map);

            if (nPartSeqLen == 1) {
                findPtmOnOneAa(nModPepsSet, flexiableMass, absPos_MassVarPtm_Map, nPartSeq, tagPosInProt-1, ms1TolAbs);
            } else {
                Map<Integer, Set<Double>> oneTimeMassGroups = new HashMap<>(nPartSeqLen);
                Map<Set<Integer>, Set<Double>> posComb_multiMassSet_Map = new HashMap<>(200);
                Map<Integer, Map<Double, VarPtm>> pos_MassVarPtm_Map = new HashMap<>(nPartSeqLen, 1);
                Map<Double, Set<Integer>> allMassAllPosesMap = new HashMap<>(nPartSeqLen * 20, 1); //Set<pos>

                inferPTM.prepareInfoNTerm(scanNum, nPartSeq, oneTimeMassGroups, posComb_multiMassSet_Map,
                        pos_MassVarPtm_Map, allMassAllPosesMap, yIdMinAbsPosMap, optStartPos == 0, fixedStartPos, tagPosInProt, protSeq);

                Set<Pair<Integer, Map<Double, Integer>>> resList = new HashSet<>(100);
                rec.x_num += inferPTM.findBestPtmMIPExtN(scanNum, env, allMassAllPosesMap, flexiableMass, tagPosInProt,
                        nPartSeq, ms1TolAbs, oneTimeMassGroups, posComb_multiMassSet_Map, posYIdMap, pos_MassVarPtm_Map, resList);
                inferPTM.getFeasibleMassPosMapN(scanNum, resList, unUsedPlMap, nPartSeq, nCutMass, N_PART,
                        unUsedExpProcessedPL, false, allMassAllPosesMap, pos_MassVarPtm_Map, nModPepsSet, tagPosInProt, yIdMinAbsPosMap, fixedStartPos);
            }
            if (yIdMinAbsPosMap.isEmpty()) {
                inferPTM.findPosssible1Ptm(scanNum, nPartSeq, absPos_MassVarPtm_Map, optStartPos, nModPepsSet, flexiableMass, ms1TolAbs);
            }
        }
        if (nModPepsSet.isEmpty()) return false;
        return true;
    }


    private void findPtmOnOneAa(TreeSet<Peptide> midModPepsSet, double deltaMass, Map<Integer, TreeMap<Double, VarPtm>> absPos_MassVarPtm_Map, String gapSeq, int absPos, double ms1TolAbs) {
        TreeMap<Double, VarPtm> massPtmMap = absPos_MassVarPtm_Map.get(absPos);
        if (massPtmMap == null) return;
        for (double ptmMass : massPtmMap.keySet()) {
            if (Math.abs(ptmMass - deltaMass) < ms1TolAbs+ms2Tolerance) {
                Peptide tmpPeptide = new Peptide(gapSeq, false, massTool);
                PosMassMap posMassMap = new PosMassMap();
                posMassMap.put(0, ptmMass);
                tmpPeptide.posVarPtmResMap.put(0, massPtmMap.get(ptmMass));
                tmpPeptide.setVarPTM(posMassMap);
                tmpPeptide.setScore(0);
                midModPepsSet.add(tmpPeptide);
            }
        }
    }

    private void updatePeptideTreeSet(Peptide newPeptide, TreeSet<Peptide> peptideTreeSet, Map<String, PeptideInfo> peptideInfoMap,
                                      String protId, String protSeq, int pepStartPos, int pepEndPos) {

        String shortProtId = protId.split(" ")[0];
        boolean added = false;
        if (peptideTreeSet.size() < candisNum) {
            peptideTreeSet.add(newPeptide);
            added = true;
        } else if (newPeptide.compareTo( peptideTreeSet.last() ) == 1) { //use the override '>' to compare because I am considering priority
            added = true;
            peptideTreeSet.pollLast();
            peptideTreeSet.add(newPeptide);
        }
        if (added) {
            char leftFlank  = pepStartPos==0 ? '-' : protSeq.charAt(pepStartPos-1);
            char rightFlank = pepEndPos==protSeq.length()-1 ? '-' : protSeq.charAt(pepEndPos+1);
            String freePepSeq = newPeptide.getFreeSeq();
            PeptideInfo peptideInfo = peptideInfoMap.get(freePepSeq);
            if (peptideInfo != null) {
                if ( ! peptideInfo.protIdSet.contains(shortProtId)) { //if this pep with prot is not recorded, add this prot
                    peptideInfo.leftFlank = leftFlank;
                    peptideInfo.rightFlank = rightFlank;
                    peptideInfo.protIdSet.add(shortProtId);
                    if (!shortProtId.startsWith("DECOY_")) {
                        peptideInfo.isDecoy = false;
                    }
                }
            } else {
                peptideInfo = new PeptideInfo(freePepSeq, shortProtId.startsWith("DECOY_"), leftFlank, rightFlank);
                peptideInfo.protIdSet.add(shortProtId);
                peptideInfoMap.put(freePepSeq, peptideInfo);
            }
        }
    }

    private void storeCleanPartPeptides(String partSeq, double cutMass, byte isNCPart, SparseVector expProcessedPL, TreeMap<Double,Double> plMap, TreeSet<Peptide> modPepsSet) {

        Peptide partPeptide = new Peptide( partSeq, false, massTool);
        double[][] ionMatrix = partPeptide.getIonMatrixNow();
        inferPTM.updateIonMatrix(ionMatrix, cutMass, isNCPart);
        Set<Integer> jRange = IntStream.rangeClosed(0, partSeq.length()-1).boxed().collect(Collectors.toSet());
        double score = massTool.buildVectorAndCalXCorr(ionMatrix, 1, expProcessedPL, partPeptide.matchedBions, partPeptide.matchedYions, jRange) ;
        partPeptide.setScore(score*(1-partPeptide.posVarPtmResMap.size()*0.01));
        partPeptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ionMatrix, ms2Tolerance));
        if (modPepsSet.size() < 2) { //max restore 2 patterns for one peptide  //todo make this avaliable to differentiate priority
            modPepsSet.add(partPeptide);
        } else if (modPepsSet.last().compareTo(partPeptide) < 0) {
            modPepsSet.pollLast();
            modPepsSet.add(partPeptide);
        }
    }

    private boolean isKR(char aa){
        return aa == 'K' || aa == 'R';
    }

    private boolean isX(char aa){
        return aa == 'X';
    }
    public class Entry {
        public Map<String, PeptideInfo> peptideInfoMapForRef = new HashMap<>();
        public Peptide topPeptide = null;
        final int scanNum;
        final String scanName;
        final int precursorCharge;
        final double precursorMass;
        final String labelling;
        public final String peptide;
        final double theoMass;
        final int isDecoy;
        public final double score;
        final double deltaLCn;
        final double deltaCn;
        final int matchedPeakNum;
        final double ionFrac;
        final double matchedHighestIntensityFrac;
        final double explainedAaFrac;
        final String otherPtmPatterns; // It has 4 decimal because it is write the the result file for checking. It is not used in scoring or other purpose.
        final String aScore;
        final String candidates;
        final String peptideSet;

        public String recString;
        List<VarPtm> varPtmList = new ArrayList<>();
        Entry(int scanNum, String scanName, int precursorCharge, double precursorMass
                ,String labelling, String peptide, double theoMass, int isDecoy
                , double score, double deltaLCn, double deltaCn, int matchedPeakNum, double ionFrac, double matchedHighestIntensityFrac
                , double explainedAaFrac, String otherPtmPatterns, String aScore, String candidates, String peptideSet) {
            this.scanNum = scanNum;
            this.scanName = scanName;
            this.precursorCharge = precursorCharge;
            this.precursorMass = precursorMass;
            this.labelling = labelling;
            this.peptide = peptide;
            this.theoMass = theoMass;
            this.isDecoy = isDecoy;
            this.score = score;
            this.deltaLCn = deltaLCn;
            this.deltaCn = deltaCn;
            this.matchedPeakNum = matchedPeakNum;
            this.ionFrac = ionFrac;
            this.matchedHighestIntensityFrac = matchedHighestIntensityFrac;
            this.explainedAaFrac = explainedAaFrac;
            this.otherPtmPatterns = otherPtmPatterns;
            this.aScore = aScore;
            this.candidates = candidates;
            this.peptideSet = peptideSet;
        }
    }
}
