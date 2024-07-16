package proteomics.PTM;

import ProteomicsLibrary.*;
import ProteomicsLibrary.MassTool;
import ProteomicsLibrary.Types.*;
import com.google.common.collect.Sets;
import gurobi.*;
import org.apache.commons.math3.util.Pair;
import proteomics.Types.*;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.dom4j.Document;
import org.dom4j.Element;
import org.dom4j.io.SAXReader;

import static proteomics.PIPI.*;


public class InferPTM {
    private static final int PoolSolutions = 100;
    private static final int MaxPtmNumInPart = 3;
    public static final byte N_PART = 0;
    public static final byte C_PART = 1;
    public static final byte ANYWHERE = 4;
    public static final byte PEPC = 3;
    public static final byte PEPN = 2;
    public static final byte PROTC = 1;
    public static final byte PROTN = 0;
    public final static DecimalFormat df3 = new DecimalFormat("0.000");
    private final MassTool massTool;
    private final Map<String, Double> elementTable;
    private final Map<Character, Double> massTable;
    private final Map<Character, Double> fixModMap;
    private Set<VarPtm> varPtmSet = new HashSet<>();
    private final double minPtmMass;
    private final double maxPtmMass;
    private final double ms2Tolerance;
    private Map<Character, List<VarPtm>> finalPtmMap = new HashMap<>(22);
    private Map<Character, Map<Byte, List<VarPtm>>> aaAllVarPtmMap = new HashMap<>(22);
    private final Set<Character> aaCharSet = new HashSet<>(Arrays.asList('A','C','D','E','F','G','H','I','K','L','M','N','O','P','Q','R','S','T','U','V','W','Y'));
    private int g_thread_num;
    private double timeLimit;
    private Set<Character> aaWithFixModSet = new HashSet<>();
    public InferPTM(MassTool massTool, Map<Character, Double> fixModMap, Map<String, String> parameterMap) throws Exception{
        this.massTool = massTool;
        elementTable = massTool.getElementTable();
        massTable = massTool.getMassTable();
        this.fixModMap = fixModMap;
        for (Character c : fixModMap.keySet()){
            if (Math.abs(fixModMap.get(c)) > 0.02) {
                aaWithFixModSet.add(c);
            }
        }
        this.minPtmMass = Double.valueOf(parameterMap.get("min_ptm_mass"));
        this.g_thread_num = Integer.valueOf(parameterMap.get("GUROBI_thread_num"));
        this.timeLimit = Double.valueOf(parameterMap.get("GUROBI_time_limit"));
        this.maxPtmMass = Double.valueOf(parameterMap.get("max_ptm_mass"));
        this.ms2Tolerance = Double.valueOf(parameterMap.get("ms2_tolerance"));

        readModFromUnimod();

        // update ptm table with the high priority mods in the parameter file.
        for (VarPtm varPtm : varPtmSet) {
            if (finalPtmMap.containsKey(varPtm.site)) {
                finalPtmMap.get(varPtm.site).add(varPtm);
            } else {
                List<VarPtm> tempList = new LinkedList<>();
                tempList.add(varPtm);
                finalPtmMap.put(varPtm.site, tempList);
            }
        }

        aaAllVarPtmMap = new HashMap<>(22); //<aachar, <position, varptm>>
        for (Character aa : finalPtmMap.keySet()){
            Map<Byte, List<VarPtm>> positionVarPtmMap = new HashMap<>(5);
            for (VarPtm varPtm : finalPtmMap.get(aa)) {
                byte position = varPtm.position;
                if (positionVarPtmMap.containsKey(position)) {
                    positionVarPtmMap.get(position).add(varPtm);
                } else {
                    List<VarPtm> varPtmList = new ArrayList<>(50);
                    varPtmList.add(varPtm);
                    positionVarPtmMap.put(position, varPtmList);
                }
            }
            aaAllVarPtmMap.put(aa, positionVarPtmMap);
        }
    }

    public void findPosssible1Ptm(int scanNum, String partSeq, Map<Integer, TreeMap<Double, VarPtm>> pos_MassVarPtm_Map, int startRefPos, TreeSet<Peptide> modPepsSet, double deltaMass, double ms1TolAbs) {

        int partSeqLen = partSeq.length();
        int absPos;
        for (int relPos = 0; relPos < partSeqLen; ++relPos) {
            absPos = relPos + startRefPos;
            TreeMap<Double, VarPtm> massVarPtmMap = pos_MassVarPtm_Map.get(absPos);
            if (massVarPtmMap == null) continue;
            for (double mass : massVarPtmMap.subMap(deltaMass-ms2Tolerance-ms1TolAbs, deltaMass+ms2Tolerance+ms1TolAbs).keySet()) {
                if (massTable.get(partSeq.charAt(relPos)) + mass < ms2Tolerance) continue;
                VarPtm varPtm = massVarPtmMap.get(mass);

                Peptide tmpPeptide = new Peptide(partSeq, false, massTool);
                PosMassMap posMassMap = new PosMassMap();
                posMassMap.put(relPos, mass);
                tmpPeptide.posVarPtmResMap.put(relPos, varPtm);
                tmpPeptide.setVarPTM(posMassMap);
                tmpPeptide.setScore(0);
                modPepsSet.add(tmpPeptide);
            }
        } // end for (int i = 0; i < partSeq.length(); ++i) {
    }
    public int findBestPtmMIPExtC(int scanNum, GRBEnv env, Map<Double, Set<Integer>> allMassAllPosesMap, double totalDeltaMass, int refPos, String partSeq,
                                   double ms1TolAbs, Map<Integer, Set<Double>> oneTimeMassGroups, final Map<Set<Integer>, Set<Double>> posComb_multiMassSet_Map,
                                   Map<Integer, Integer> posYIdMap, Map<Integer, Map<Double, VarPtm>> absPos_MassVarPtm_Map, Set<Pair<Integer, Map<Double, Integer>>> resList) {
        int varNum = 0;
        Map<Integer, List<Integer>> yIdAllPosesMap = new HashMap<>(posYIdMap.values().size());
        for (int pos : posYIdMap.keySet()) {
            int yId = posYIdMap.get(pos);
            List<Integer> allPoses = yIdAllPosesMap.get(yId);
            if ( allPoses != null) {
                allPoses.add(pos);
            } else {
                allPoses = new ArrayList<>();
                allPoses.add(pos);
                yIdAllPosesMap.put(yId, allPoses);
            }
        }
        try {
            GRBModel model = new GRBModel(env);
            //// Constraints
            GRBLinExpr totalNumsOnPepConstr = new GRBLinExpr();
            GRBLinExpr totalFlexiableMass = new GRBLinExpr();
            GRBLinExpr totalFlexiableMass1 = new GRBLinExpr();
            GRBVar t = model.addVar(0, ms1TolAbs, 1, GRB.CONTINUOUS, "t");
            totalFlexiableMass.addTerm(-1,t);
            totalFlexiableMass1.addTerm(1,t);
            // x y variables
            Map<Integer, GRBVar> yVarMap = new HashMap<>( posYIdMap.values().size() ); // <yId, var>
            for (int yId : posYIdMap.values()) {
                GRBVar yVar = model.addVar(0, 1, 0.0, GRB.BINARY, "y"+yId);
                yVarMap.put(yId, yVar);

                //constraints
                double aaMass = 0;
                for (int absPos : yIdAllPosesMap.get(yId)) {
                    aaMass += massTool.getMassTable().get(partSeq.charAt(absPos-refPos));
                }
                totalFlexiableMass.addTerm(aaMass, yVarMap.get(yId));
                totalFlexiableMass1.addTerm(aaMass, yVarMap.get(yId));

            }
            List<Integer> yIdList = new ArrayList<>(yVarMap.keySet());
            yIdList.sort(Comparator.naturalOrder());
            Map<Double, GRBVar> xVarMap = new HashMap<>(allMassAllPosesMap.size());
            for (double mass : allMassAllPosesMap.keySet()) {
                Set<Integer> allAbsPoses = allMassAllPosesMap.get(mass);
                GRBVar xVar = model.addVar(0, allAbsPoses.size(), 0.1, GRB.INTEGER, "x_"+mass);
                xVarMap.put(mass, xVar);

                //constraints
                totalFlexiableMass.addTerm(mass, xVar); // += m_i * x_i
                totalFlexiableMass1.addTerm(mass, xVar); // += m_i * x_i

                totalNumsOnPepConstr.addTerm(1,xVar); // + 1 * x_i

                //constraints
                GRBLinExpr massOccurence = new GRBLinExpr();
                massOccurence.addTerm(1, xVar);
                int fixPosNum = 0;

                for (int absPos : allAbsPoses) {
                    if ( ! posYIdMap.containsKey(absPos)) {
                        fixPosNum++;
                        continue;
                    }

                    int yId = posYIdMap.get(absPos);
                    int position = absPos_MassVarPtm_Map.get(absPos).get(mass).position;

                    if (position == PEPC
                            && absPos-refPos != partSeq.length()-1
                            && yId != yIdList.get(yIdList.size()-1)) {
                        massOccurence.addTerm(-1, yVarMap.get(yId)); // thisYId
                        massOccurence.addTerm(1, yVarMap.get(yId+1));// nextYId
                    } else {
                        massOccurence.addTerm(-1, yVarMap.get(yId));
                    }
                }
                if (massOccurence.size() > 1) { //if it is just x <= 1, no any y, then no need to add this constraint
                    model.addConstr(massOccurence, GRB.LESS_EQUAL, fixPosNum, "massOccurence_"+mass); // xi <= fixPosNum + y1+y2+...
                }
            }

            //// add constraints
            model.addConstr(totalFlexiableMass1, GRB.GREATER_EQUAL, totalDeltaMass , "totalFlexiableMassGe");
            model.addConstr(totalFlexiableMass, GRB.LESS_EQUAL, totalDeltaMass, "totalFlexiableMassLe");
            model.addConstr(totalNumsOnPepConstr, GRB.LESS_EQUAL, Math.min(partSeq.length(), MaxPtmNumInPart), "totalNumsOnPepConstr"); // this value should not exceed the number of aa in the partSeq

            // constraints, Sum(oneMass on certain aa) < 1 or y
            for (int absPos : oneTimeMassGroups.keySet()) {
                GRBLinExpr sumX_leq_1orY = new GRBLinExpr();
//                int absPos = refPos + pos;
                for (double mass : oneTimeMassGroups.get(absPos)) {
                    sumX_leq_1orY.addTerm(1, xVarMap.get(mass));
                }
                if (posYIdMap.containsKey(absPos)){
                    sumX_leq_1orY.addTerm(-1, yVarMap.get(posYIdMap.get(absPos)));
                    model.addConstr(sumX_leq_1orY, GRB.LESS_EQUAL, 0, "sumX_leq_Y" + absPos);
                } else {
                    model.addConstr(sumX_leq_1orY, GRB.LESS_EQUAL, 1, "sumX_leq_1" + absPos);
                }
            }

            int dummyI = 0;
            for (Set<Integer> posComb : posComb_multiMassSet_Map.keySet()) {
                GRBLinExpr multiTimeMassConstr = new GRBLinExpr();
                int fixPosNum = 0;
                for (int absPos : posComb) {
                    if (posYIdMap.containsKey(absPos)) {
                        multiTimeMassConstr.addTerm(-1, yVarMap.get(posYIdMap.get(absPos)));
                    } else {
                        fixPosNum++;
                    }
                }
                Set<Double> massSet = posComb_multiMassSet_Map.get(posComb);
                for (double mass : massSet) {
                    multiTimeMassConstr.addTerm(1, xVarMap.get(mass));
                }
                model.addConstr(multiTimeMassConstr, GRB.LESS_EQUAL, fixPosNum, "multiTimeMassConstr_"+fixPosNum+"_"+dummyI);
                dummyI++;
            }

            if (yVarMap.size() >= 2) {
                for (int i = 0; i < yIdList.size()-1; i++) {
                    int thisYId = yIdList.get(i);
                    int nextYId = yIdList.get(i+1);
                    GRBLinExpr yiyj = new GRBLinExpr();
                    yiyj.addTerm(1, yVarMap.get(thisYId));
                    yiyj.addTerm(-1, yVarMap.get(nextYId));
                    model.addConstr(yiyj, GRB.GREATER_EQUAL, 0 , "yiyj_"+i);
                }
            }

            varNum = xVarMap.size() + yVarMap.size();
            //obj function
//            model.setObjective(totalNumsOnPepConstr, GRB.MINIMIZE); // with this, the solver will find those with smaller number of PTM first then more numbers
//            GRBLinExpr tConstr = new GRBLinExpr();
//            tConstr.addTerm(1, t);
//            model.setObjective(tConstr, GRB.MINIMIZE); // with this, the solver will find those with smaller number of PTM first then more numbers
            model.set(GRB.IntAttr.ModelSense, GRB.MINIMIZE);
            if (debugScanNum.contains(scanNum) && partSeq.contentEquals("KFGVLSDNFK")){
                int a = 1;
            }
            model.set(GRB.IntParam.MIPFocus, 1); // 2 seems better than 1 but dont know why
            model.set(GRB.IntParam.PoolSearchMode, 2 );  //0 for only one sol, 1 for possible more but not guaranteed = poolSolutions, 2 for guaranteed but long time
            model.set(GRB.IntParam.PoolSolutions, PoolSolutions);
            model.set(GRB.DoubleParam.TimeLimit, timeLimit); // second
            model.set(GRB.IntParam.ConcurrentMIP, 2); // second
            model.set(GRB.IntParam.Threads, g_thread_num); // second

            model.optimize();
            switch (model.get(GRB.IntAttr.Status)) {
                case GRB.OPTIMAL:
                    break;
                case GRB.TIME_LIMIT:
                    break;
                default:
                    model.dispose();
                    return varNum; // if it is disposed here (i.e. infeasible, unbounded...), dont collect solutions, just return
            }
            int solCount = model.get(GRB.IntAttr.SolCount);
            for (int solNum = 0; solNum < solCount; solNum++) {
                model.set(GRB.IntParam.SolutionNumber, solNum);

                int maxYId = -99;
                for (int yId : yVarMap.keySet()) {
                    GRBVar yVar = yVarMap.get(yId);
                    int time = (int) Math.round(yVar.get(GRB.DoubleAttr.Xn));
                    if (time > 0) {
                        if (yId > maxYId) {
                            maxYId = yId;
                        }
                    }
                }
                Map<Double, Integer> massTimeMap = new HashMap<>();
                for (double mass : xVarMap.keySet()){
                    GRBVar xVar = xVarMap.get(mass);
                    int time = (int) Math.round(xVar.get(GRB.DoubleAttr.Xn));
                    if (time > 0) {
                        massTimeMap.put(mass, time);
                    }
                }
                resList.add(new Pair<>(maxYId, massTimeMap));
            }
            model.dispose();
        } catch (GRBException e) {
            System.out.println(scanNum + "," + partSeq+ " , Error code: " + e.getErrorCode() + ". " + e.getMessage());
        }
        return varNum;
    }

    public int findBestPtmMIPExtN(int scanNum, GRBEnv env, Map<Double, Set<Integer>> allMassAllPosesMap, double totalDeltaMass, int refPos, String partSeq,
                                   double ms1TolAbs, Map<Integer, Set<Double>> oneTimeMassGroups, final Map<Set<Integer>, Set<Double>> posComb_multiMassSet_Map,
                                   Map<Integer, Integer> posYIdMap, Map<Integer, Map<Double, VarPtm>> absPos_MassVarPtm_Map, Set<Pair<Integer, Map<Double, Integer>>> resList) {

        int varNum = 0;
        Map<Integer, List<Integer>> yIdAllPosesMap = new HashMap<>();
        for (int pos : posYIdMap.keySet()) {
            int yId = posYIdMap.get(pos);
            List<Integer> allPoses = yIdAllPosesMap.get(yId);
            if ( allPoses != null) {
                allPoses.add(pos);
            } else {
                allPoses = new ArrayList<>();
                allPoses.add(pos);
                yIdAllPosesMap.put(yId, allPoses);
            }
        }

        try {
            GRBModel model = new GRBModel(env);
            //// Constraints
            GRBLinExpr totalNumsOnPepConstr = new GRBLinExpr();
            GRBLinExpr totalFlexiableMass = new GRBLinExpr();
            GRBLinExpr totalFlexiableMass1 = new GRBLinExpr();
            GRBVar t = model.addVar(0, ms1TolAbs, 1, GRB.CONTINUOUS, "t");
            totalFlexiableMass.addTerm(-1,t);
            totalFlexiableMass1.addTerm(1,t);
            // x y variables
            Map<Integer, GRBVar> yVarMap = new HashMap<>( posYIdMap.values().size() ); // <yId, var>
            for (int yId : posYIdMap.values()) {
                GRBVar yVar = model.addVar(0, 1, 0.0, GRB.BINARY, "y"+yId);
                yVarMap.put(yId, yVar);

                //constraints
                double aaMass = 0;
                for (int absPos : yIdAllPosesMap.get(yId)) {
                    aaMass += massTool.getMassTable().get(partSeq.charAt(absPos + partSeq.length()-refPos));
                }
                totalFlexiableMass.addTerm(aaMass, yVarMap.get(yId));
                totalFlexiableMass1.addTerm(aaMass, yVarMap.get(yId));

            }
            List<Integer> yIdList = new ArrayList<>(yVarMap.keySet());
            yIdList.sort(Comparator.naturalOrder());
            Map<Double, GRBVar> xVarMap = new HashMap<>(allMassAllPosesMap.size());
            for (double mass : allMassAllPosesMap.keySet()) {
                Set<Integer> allAbsPoses = allMassAllPosesMap.get(mass);
                GRBVar xVar = model.addVar(0, allAbsPoses.size(), 0.1, GRB.INTEGER, "x_"+mass);
                xVarMap.put(mass, xVar);

                //constraints
                totalFlexiableMass.addTerm(mass, xVar); // += m_i * x_i
                totalFlexiableMass1.addTerm(mass, xVar); // += m_i * x_i

                totalNumsOnPepConstr.addTerm(1,xVar); // + 1 * x_i

                //constraints
                GRBLinExpr massOccurence = new GRBLinExpr();
                massOccurence.addTerm(1, xVar);
                int fixPosNum = 0;

                for (int absPos : allAbsPoses) {
                    if ( ! posYIdMap.containsKey(absPos)) {
                        fixPosNum++;
                        continue;
                    }

                    int yId = posYIdMap.get(absPos);
                    int position = absPos_MassVarPtm_Map.get(absPos).get(mass).position;

                    if (position == PEPN
                            && absPos+partSeq.length() != refPos
                            && yId != yIdList.get(yIdList.size()-1)) { // NC the same
                        massOccurence.addTerm(-1, yVarMap.get(yId)); // thisYId
                        massOccurence.addTerm(1, yVarMap.get(yId+1));// nextYId // NC the same
                    } else {
                        massOccurence.addTerm(-1, yVarMap.get(yId));
                    }
                }
                if (massOccurence.size() > 1) { //if it is just x <= 1, no any y, then no need to add this constraint
                    model.addConstr(massOccurence, GRB.LESS_EQUAL, fixPosNum, "massOccurence_"+mass); // xi <= fixPosNum + y1+y2+...
                }
            }

            //// add constraints
            model.addConstr(totalFlexiableMass1, GRB.GREATER_EQUAL, totalDeltaMass , "totalFlexiableMassGe");
            model.addConstr(totalFlexiableMass, GRB.LESS_EQUAL, totalDeltaMass, "totalFlexiableMassLe");
            model.addConstr(totalNumsOnPepConstr, GRB.LESS_EQUAL, Math.min(partSeq.length(), MaxPtmNumInPart), "totalNumsOnPepConstr"); // this value should not exceed the number of aa in the partSeq

            // constraints, Sum(oneMass on certain aa) < 1 or y
            for (int absPos : oneTimeMassGroups.keySet()) {
                GRBLinExpr sumX_leq_1orY = new GRBLinExpr();
                for (double mass : oneTimeMassGroups.get(absPos)) {
                    sumX_leq_1orY.addTerm(1, xVarMap.get(mass));
                }
                if (posYIdMap.containsKey(absPos)){
                    sumX_leq_1orY.addTerm(-1, yVarMap.get(posYIdMap.get(absPos)));
                    model.addConstr(sumX_leq_1orY, GRB.LESS_EQUAL, 0, "sumX_leq_Y" + absPos);
                } else {
                    model.addConstr(sumX_leq_1orY, GRB.LESS_EQUAL, 1, "sumX_leq_1" + absPos);
                }
            }

            int dummyI = 0;
            for (Set<Integer> posComb : posComb_multiMassSet_Map.keySet()) {
                GRBLinExpr multiTimeMassConstr = new GRBLinExpr();
                int fixPosNum = 0;
                for (int absPos : posComb) {
                    if (posYIdMap.containsKey(absPos)) {
                        multiTimeMassConstr.addTerm(-1, yVarMap.get(posYIdMap.get(absPos)));
                    } else {
                        fixPosNum++;
                    }
                }
                Set<Double> massSet = posComb_multiMassSet_Map.get(posComb);
                for (double mass : massSet) {
                    multiTimeMassConstr.addTerm(1, xVarMap.get(mass));
                }
                model.addConstr(multiTimeMassConstr, GRB.LESS_EQUAL, fixPosNum, "multiTimeMassConstr_"+fixPosNum+"_"+dummyI);
                dummyI++;
            }

            if (yVarMap.size() >= 2) {
                for (int i = 0; i < yIdList.size()-1; i++) {
                    int thisYId = yIdList.get(i);
                    int nextYId = yIdList.get(i+1);
                    GRBLinExpr yiyj = new GRBLinExpr();
                    yiyj.addTerm(1, yVarMap.get(thisYId));
                    yiyj.addTerm(-1, yVarMap.get(nextYId));
                    model.addConstr(yiyj, GRB.GREATER_EQUAL, 0 , "yiyj_"+i);
                }
            }

            varNum = xVarMap.size() + yVarMap.size();
            //obj function
//            GRBLinExpr tConstr = new GRBLinExpr();
//            tConstr.addTerm(1, t);
//            model.setObjective(tConstr, GRB.MINIMIZE); // with this, the solver will find those with smaller number of PTM first then more numbers
            model.set(GRB.IntAttr.ModelSense, GRB.MINIMIZE);

            if (debugScanNum.contains(scanNum) && partSeq.contentEquals("KFGVLSDNFK")){
                int a = 1;
            }
            model.set(GRB.IntParam.MIPFocus, 1); // 2 seems better than 1 but dont know why
            model.set(GRB.IntParam.PoolSearchMode, 2 );  //0 for only one sol, 1 for possible more but not guaranteed = poolSolutions, 2 for guaranteed but long time
            model.set(GRB.IntParam.PoolSolutions, PoolSolutions);
            model.set(GRB.DoubleParam.TimeLimit, timeLimit); // second
            model.set(GRB.IntParam.ConcurrentMIP, 2); // second
            model.set(GRB.IntParam.Threads, g_thread_num); // second

            model.optimize();
            switch (model.get(GRB.IntAttr.Status)) {
                case GRB.OPTIMAL:
                    break;
                case GRB.TIME_LIMIT:
                    break;
                default:
                    model.dispose();
                    return varNum;
            }
            int solCount = model.get(GRB.IntAttr.SolCount);
            for (int solNum = 0; solNum < solCount; solNum++) {
                model.set(GRB.IntParam.SolutionNumber, solNum);

                int maxYId = -99;
                for (int yId : yVarMap.keySet()) {
                    GRBVar yVar = yVarMap.get(yId);
                    int time = (int) Math.round(yVar.get(GRB.DoubleAttr.Xn));
                    if (time > 0) {
                        if (yId > maxYId) {
                            maxYId = yId;
                        }
                    }
                }
                Map<Double, Integer> massTimeMap = new HashMap<>();
                for (double mass : xVarMap.keySet()){
                    GRBVar xVar = xVarMap.get(mass);
                    int time = (int) Math.round(xVar.get(GRB.DoubleAttr.Xn));
                    if (time > 0) {
                        massTimeMap.put(mass, time);
                    }
                }
                resList.add(new Pair<>(maxYId, massTimeMap));
            }
            model.dispose();
        } catch (GRBException e) {
            System.out.println(scanNum + "," + partSeq+ " , Error code: " + e.getErrorCode() + ". " + e.getMessage());
        }
        return varNum;
    }

    public void updateIonMatrix(double[][] ionMatrix, double cutMass, byte ncPart){
        if (ncPart == N_PART) { // is n part seq
            for (int i = 0; i < ionMatrix[1].length; i++) {
                ionMatrix[1][i] += cutMass;
            }
        } else {            //is c part seq
            for (int i = 0; i < ionMatrix[0].length; i++) {
                ionMatrix[0][i] += cutMass;
            }
        }
    }

    public void getFeasibleMassPosMapC(int scanNum, Set<Pair<Integer, Map<Double, Integer>>> massTimeResList, TreeMap<Double, Double> plMap, String fullPartSeq,
                                       double cutMass, byte ncPart, SparseVector expProcessedPL, boolean isDecoy, Map<Double, Set<Integer>> allMassAllPosesMap,
                                       Map<Integer, Map<Double, VarPtm>> absPos_MassVarPtm_Map, TreeSet<Peptide> cModPepsSet, int startRefPos, Map<Integer, Integer> yIdMaxAbsPosMap,
                                       int optStartPos) {

        Map<Integer, Double> singleAbsPosMassMap;
        Map<Double, Integer> massTimeMap;
        for (Pair<Integer, Map<Double, Integer>> y_x_Res : massTimeResList) {
            int maxYId = y_x_Res.getFirst();
            int seqEndPos;

            if (maxYId == -99) {
                seqEndPos = optStartPos;
            } else {
                seqEndPos = yIdMaxAbsPosMap.get(maxYId)+1;
            }

            String partSeq = fullPartSeq.substring(0, seqEndPos-startRefPos);

            massTimeMap = y_x_Res.getSecond();
            singleAbsPosMassMap = new HashMap<>();
            List<Double> massMaxMultiTime = new ArrayList<>(); // the mass order in this list is fixed. should be reused to match there poses.
            for (double mass : massTimeMap.keySet()){
                Set<Integer> allTheoAbsPoses = allMassAllPosesMap.get(mass);
                if (allTheoAbsPoses.size() == 1) { // max poses num is 1. First use them to occupy some AAs
                    singleAbsPosMassMap.put(Collections.min(allTheoAbsPoses), mass); // not finding minimum just get the element
                } else {
                    for (int i = 0; i < massTimeMap.get(mass); i++) {  // if in res, one mass has xi = 2, then add it twice for correct cartesianProduct later
                        massMaxMultiTime.add(mass);
                    }                }
            } // finish checking all res mass whose MaxPosSize==1


            List<Set<Integer>> updatedFeasiblePosesList = new ArrayList<>(massMaxMultiTime.size());
            for (double mass : massMaxMultiTime) { // for mass whose MaxPosSize > 1
                Set<Integer> feasibleAbsPoses = new HashSet<>();
                for (int absPos : allMassAllPosesMap.get(mass)) {
                    if (absPos < seqEndPos && !singleAbsPosMassMap.containsKey(absPos)) {
                        feasibleAbsPoses.add(absPos);
                    }
                }
                updatedFeasiblePosesList.add(feasibleAbsPoses);
            }

            Set<PosMassMap> triedPtmPattern = new HashSet<>(200);
            // then the updatedFeasiblePosesList is empty, suprisingly it still works to find the only-single-time PTM pattern.
            for (List<Integer> cartesianList : Sets.cartesianProduct(updatedFeasiblePosesList)){
                Set<Integer> cartesianSet = new HashSet<>(cartesianList);
                if (cartesianSet.size() < cartesianList.size()) continue; // if the cartesian list has repeated elements, i.e., more than one ptm on one pos, skip

                Peptide tmpPeptide = new Peptide(partSeq, isDecoy, massTool);
                PosMassMap posMassMap = new PosMassMap();

                for (int absPos : singleAbsPosMassMap.keySet()) {
                    int relPos = absPos - startRefPos;
                    posMassMap.put(relPos, singleAbsPosMassMap.get(absPos)); // set every the single Pos Ptms
                    tmpPeptide.posVarPtmResMap.put(relPos, absPos_MassVarPtm_Map.get(absPos).get(singleAbsPosMassMap.get(absPos)));
                }
                for (int i = 0; i < massMaxMultiTime.size(); i++) {
                    double mass = massMaxMultiTime.get(i);
                    int absPos = cartesianList.get(i);
                    int relPos = absPos - startRefPos;
                    posMassMap.put(relPos, mass); // set every the single Pos Ptms
                    tmpPeptide.posVarPtmResMap.put(relPos, absPos_MassVarPtm_Map.get(absPos).get(mass));
                }

                if (triedPtmPattern.contains(posMassMap)) {
                    continue;
                }else {
                    triedPtmPattern.add(posMassMap);
                }
                double[][] ionMatrix = tmpPeptide.getIonMatrix();
                if ( ! posMassMap.isEmpty()) { // if it is empty, the ptm must be just proton and ignored, so we just use the original ionMatrix to calculate score
                    tmpPeptide.setVarPTM(posMassMap);
                    ionMatrix = tmpPeptide.getIonMatrixNow();
                    updateIonMatrix(ionMatrix, cutMass, ncPart);
                }

                Set<Integer> jRange = IntStream.rangeClosed(0, partSeq.length()-1).boxed().collect(Collectors.toSet());
                double score = massTool.buildVectorAndCalXCorr(ionMatrix, 1, expProcessedPL, tmpPeptide.matchedBions, tmpPeptide.matchedYions, jRange) ;
                tmpPeptide.setScore(score*(1-tmpPeptide.posVarPtmResMap.size()*0.05));
                tmpPeptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, tmpPeptide.getIonMatrix(), ms2Tolerance));
                if (cModPepsSet.size() < 2) { //max restore 2 patterns for one peptide  //todo make this avaliable to differentiate priority
                    cModPepsSet.add(tmpPeptide);
                } else if (cModPepsSet.last().compareTo(tmpPeptide) < 0) {
                    cModPepsSet.pollLast();
                    cModPepsSet.add(tmpPeptide);
                }
            }
        }
        int a = 1;
    }

    public void getFeasibleMassPosMapN(int scanNum, Set<Pair<Integer, Map<Double, Integer>>> massTimeResList, TreeMap<Double, Double> plMap, String fullPartSeq,
                                       double cutMass, byte ncPart, SparseVector expProcessedPL, boolean isDecoy, Map<Double, Set<Integer>> allMassAllPosesMap,
                                       Map<Integer, Map<Double, VarPtm>> absPos_MassVarPtm_Map, TreeSet<Peptide> nModPepsSet, int startRefPos, Map<Integer, Integer> yIdMinAbsPosMap,
                                       int maxNPos) {

        int maxYId;
        int seqStartPos;
        String partSeq;
        Map<Double, Integer> massTimeMap;
        Map<Integer, Double> singleAbsPosMassMap;
        List<Double> massMaxMultiTime;
        List<Set<Integer>> updatedFeasiblePosesList;
        for (Pair<Integer, Map<Double, Integer>> y_x_Res : massTimeResList) {
            maxYId = y_x_Res.getFirst();
            if (maxYId == -99) {
                seqStartPos = maxNPos;
            } else {
                seqStartPos = yIdMinAbsPosMap.get(maxYId);
            }

            partSeq = fullPartSeq.substring(seqStartPos-startRefPos + fullPartSeq.length());

            massTimeMap = y_x_Res.getSecond();
            singleAbsPosMassMap = new HashMap<>();
            massMaxMultiTime = new ArrayList<>(); // the mass order in this list is fixed. should be reused to match there poses.
            Set<Integer> allTheoAbsPoses;
            for (double mass : massTimeMap.keySet()){
                allTheoAbsPoses = allMassAllPosesMap.get(mass);
                if (allTheoAbsPoses.size() == 1) { // max poses num is 1. First use them to occupy some AAs
                    singleAbsPosMassMap.put(Collections.min(allTheoAbsPoses), mass);
                } else {
                    for (int i = 0; i < massTimeMap.get(mass); i++) {  // if in res, one mass has xi = 2, then add it twice for correct cartesianProduct later
                        massMaxMultiTime.add(mass);
                    }                }
            } // finish checking all res mass whose MaxPosSize==1

            updatedFeasiblePosesList = new ArrayList<>(massMaxMultiTime.size());
            for (double mass : massMaxMultiTime) { // for mass whose MaxPosSize > 1
                Set<Integer> feasibleAbsPoses = new HashSet<>();
                for (int absPos : allMassAllPosesMap.get(mass)) {
                    if (absPos >= seqStartPos && !singleAbsPosMassMap.containsKey(absPos)) {
                        feasibleAbsPoses.add(absPos);
                    }
                }
                updatedFeasiblePosesList.add(feasibleAbsPoses);
            }

            Set<PosMassMap> triedPtmPattern = new HashSet<>(200);
            // then the updatedFeasiblePosesList is empty, suprisingly it still works to find the only-single-time PTM pattern.
            for (List<Integer> cartesianList : Sets.cartesianProduct(updatedFeasiblePosesList)){
                Set<Integer> cartesianSet = new HashSet<>(cartesianList);
                if (cartesianSet.size() < cartesianList.size()) continue; // if the cartesian list has repeated elements, i.e., more than one ptm on one pos, skip

                Peptide tmpPeptide = new Peptide(partSeq, isDecoy, massTool);
                PosMassMap posMassMap = new PosMassMap();

                for (int absPos : singleAbsPosMassMap.keySet()) {
                    int relPos = absPos - startRefPos + partSeq.length();
                    posMassMap.put(relPos, singleAbsPosMassMap.get(absPos)); // set every the single Pos Ptms
                    tmpPeptide.posVarPtmResMap.put(relPos, absPos_MassVarPtm_Map.get(absPos).get(singleAbsPosMassMap.get(absPos)));
                }
                for (int i = 0; i < massMaxMultiTime.size(); i++) {
                    double mass = massMaxMultiTime.get(i);
                    int absPos = cartesianList.get(i);
                    int relPos = absPos - startRefPos + partSeq.length();
                    posMassMap.put(relPos, mass); // set every the single Pos Ptms
                    tmpPeptide.posVarPtmResMap.put(relPos, absPos_MassVarPtm_Map.get(absPos).get(mass));
                }

                if (triedPtmPattern.contains(posMassMap)) {
                    continue;
                }else {
                    triedPtmPattern.add(posMassMap);
                }
                double[][] ionMatrix = tmpPeptide.getIonMatrix();
                if ( ! posMassMap.isEmpty()) { // if it is empty, the ptm must be just proton and ignored, so we just use the original ionMatrix to calculate score
                    tmpPeptide.setVarPTM(posMassMap);
                    ionMatrix = tmpPeptide.getIonMatrixNow();
                    updateIonMatrix(ionMatrix, cutMass, ncPart);
                }

                Set<Integer> jRange = IntStream.rangeClosed(0, partSeq.length()-1).boxed().collect(Collectors.toSet());
                double score = massTool.buildVectorAndCalXCorr(ionMatrix, 1, expProcessedPL, tmpPeptide.matchedBions, tmpPeptide.matchedYions, jRange) ;
//                if (score > 0) {
                tmpPeptide.setScore(score*(1-tmpPeptide.posVarPtmResMap.size()*0.05));
                tmpPeptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, tmpPeptide.getIonMatrix(), ms2Tolerance));
                if (nModPepsSet.size() < 2) { //max restore 5 patterns for one peptide
                    nModPepsSet.add(tmpPeptide);
                } else if (nModPepsSet.last().compareTo(tmpPeptide) < 0) {
                    nModPepsSet.pollLast();
                    nModPepsSet.add(tmpPeptide);
                }
//                }
            }
        }
        int a = 1;
    }


    public double getMinPtmMass() {
        return minPtmMass;
    }

    public double getMaxPtmMass() {
        return maxPtmMass;
    }

    private void readModFromUnimod() throws Exception {
        SAXReader reader = new SAXReader();
        InputStream inputStream = getClass().getClassLoader().getResourceAsStream("unimod.xml"); // PTMs from Unimod except for AA substitutions, isotopic labellings
        Document document = reader.read(inputStream);
        Element rootElement = document.getRootElement();
        Iterator<Element> rootIter = rootElement.elementIterator();

        while (rootIter.hasNext()) {
            Element rootElem = rootIter.next();
            if (!rootElem.getName().contentEquals("modifications")) continue;

            Iterator<Element> modIter = rootElem.elementIterator();

            while (modIter.hasNext()) {
                Element modElem = modIter.next();

                String name = modElem.attributeValue("title");
                double mass = Double.valueOf(modElem.element("delta").attributeValue("mono_mass"));
                if (mass < minPtmMass || mass > maxPtmMass) continue;
                for (Element spec : modElem.elements("specificity")) {
                    String classification = spec.attributeValue("classification");
                    if ( classification.contains("glycos") || classification.contains("Other")) {
                        continue;
                    }
                    String siteStr = spec.attributeValue("site");
                    String positionStr = spec.attributeValue("position");
                    if (classification.contentEquals("Isotopic label") && !(name.contentEquals("Propionyl") && siteStr.contentEquals("K")) && !(name.contentEquals("Succinyl") && siteStr.contentEquals("K"))) { // only for synthetic ptm data, because the authors uses them
                        continue;
                    }
                    byte position = 0;
                    switch (positionStr) {
                        case "Protein N-term":
                            position = PROTN;
                            break;
                        case "Protein C-term":
                            position = PROTC;
                            break;
                        case "Any N-term":
                            position = PEPN;
                            break;
                        case "Any C-term":
                            position = PEPC;
                            break;
                        case "Anywhere":
                            position = ANYWHERE;
                            break;
                    }
                    if (siteStr.contentEquals("N-term") || siteStr.contentEquals("C-term")) {
                        for (char site : aaCharSet) {
                            if (aaWithFixModSet.contains(site) && siteStr.contentEquals("C-term")) {
                                continue; // if aa is C, just ignore the mod that are not at N term
                            }

                            VarPtm temp = new VarPtm(mass, site, position, name, classification, 0);
                            if (finalPtmMap.containsKey(site)) {
                                finalPtmMap.get(site).add(temp);
                            } else {
                                List<VarPtm> varPtmSet = new LinkedList<>();
                                varPtmSet.add(temp);
                                finalPtmMap.put(site, varPtmSet);
                            }
                        }
                    } else {
                        char site = siteStr.charAt(0);
                        if (aaWithFixModSet.contains(site) && (position == 1 || position == 3 || position == 4)) {
                            continue;  // if aa is C, just ignore the mod that are not at N term
                        }
                        VarPtm temp = new VarPtm(mass, site, position, name, classification, 0);
                        if (finalPtmMap.containsKey(site)) {
                            finalPtmMap.get(site).add(temp);
                        } else {
                            List<VarPtm> varPtmSet = new LinkedList<>();
                            varPtmSet.add(temp);
                            finalPtmMap.put(site, varPtmSet);
                        }
                    }
                }
            }
        }
    }

    public void prepareInfoCTerm(int scanNum, String partSeq,
                                 Map<Integer, TreeMap<Double, VarPtm>> pos_MassVarPtm_Map,
                                 final Map<Integer, Integer> yIdMaxAbsPosMap,
                                 final boolean couldBeProtC,
                                 final int optStartPos,
                                 final int startRefPos,
                                 Map<Integer, List<Byte>> absPos_ptmPositions_Map) {

        int partSeqLen = partSeq.length();
        char aa;
        int absPos;
        for (int relPos = 0; relPos < partSeqLen; ++relPos) {
            aa = partSeq.charAt(relPos);
            absPos = relPos + startRefPos;
            List<Byte> positionsToTry = new ArrayList<>(3);
            absPos_ptmPositions_Map.put(absPos, positionsToTry);
            if (aaWithFixModSet.contains(aa)) continue;

            positionsToTry.add(ANYWHERE);
            if (yIdMaxAbsPosMap.containsValue(absPos)) {
                positionsToTry.add(PEPC);
            }
            if (absPos == optStartPos-1) {
                if ((aa == 'K' || aa == 'R') || !cTermSpecific) {
                    positionsToTry.add(PEPC);
                }
            }
            if (couldBeProtC && relPos == partSeqLen-1) {
                positionsToTry.add(PROTC);
            }

            if (aaAllVarPtmMap.containsKey(aa)) {
                Map<Byte, List<VarPtm>> allVarPtmMap = aaAllVarPtmMap.get(aa);
                Map<String, VarPtm> dstMap = new HashMap<>(128);
                for (Byte position : positionsToTry) {
                    if ( ! allVarPtmMap.containsKey(position)) continue;
                    for (VarPtm varPtm : allVarPtmMap.get(position)) {
                        double mass = varPtm.mass;
                        if (massTable.get(partSeq.charAt(relPos)) + mass < ms2Tolerance) continue;
                        String varPtmStr = varPtm.getStr();

                        VarPtm oldVarPtm = dstMap.get(varPtmStr);
                        if (oldVarPtm != null) {
                            if (varPtm.priority > oldVarPtm.priority){
                                dstMap.put(varPtmStr, varPtm);
                            }
                        } else {
                            dstMap.put(varPtmStr, varPtm);
                        }
                    }
                }
                if (!dstMap.isEmpty()) {
                    TreeMap<Double, VarPtm> massVarPtmMap = new TreeMap<>();
                    for (VarPtm varPtm : dstMap.values()){
                        massVarPtmMap.put(varPtm.mass, varPtm);
                    }
                    pos_MassVarPtm_Map.put(absPos, massVarPtmMap);
                }
            }
        } // end for (int i = 0; i < partSeq.length(); ++i) {
    }
    public void prepareInfoCTerm(int scanNum, String partSeq,
                                 Map<Integer, Set<Double>> oneTimeMassSetsMap,
                                 Map<Set<Integer>, Set<Double>> posComb_multiMassSet_Map,
                                 Map<Integer, Map<Double, VarPtm>> pos_MassVarPtm_Map,
                                 Map<Double, Set<Integer>> allMassAllPosesMap,
                                 final Map<Integer, Integer> yIdMaxAbsPosMap,
                                 final boolean couldBeProtC,
                                 final int optStartPos,
                                 final int startRefPos) {

        int partSeqLen = partSeq.length();
        char aa;
        int absPos;
        for (int relPos = 0; relPos < partSeqLen; ++relPos) {
            aa = partSeq.charAt(relPos);
            absPos = relPos + startRefPos;
            if (aaWithFixModSet.contains(aa)) continue;

            List<Byte> positionsToTry = new ArrayList<>(3);
            positionsToTry.add(ANYWHERE);
            if (yIdMaxAbsPosMap.containsValue(absPos)) {
                positionsToTry.add(PEPC);
            }
            if (absPos == optStartPos-1) {
                if ((aa == 'K' || aa == 'R') || !cTermSpecific) {
                    positionsToTry.add(PEPC);
                }
            }
            if (couldBeProtC && relPos == partSeqLen-1) {
                positionsToTry.add(PROTC);
            }

            if (aaAllVarPtmMap.containsKey(aa)) {
                Map<Byte, List<VarPtm>> allVarPtmMap = aaAllVarPtmMap.get(aa);
                Map<String, VarPtm> dstMap = new HashMap<>(128);
                for (Byte position : positionsToTry) {
                    if ( ! allVarPtmMap.containsKey(position)) continue;
                    for (VarPtm varPtm : allVarPtmMap.get(position)) {
                        double mass = varPtm.mass;
                        if (massTable.get(partSeq.charAt(relPos)) + mass < ms2Tolerance) continue;
                        String varPtmStr = varPtm.getStr();

                        VarPtm oldVarPtm = dstMap.get(varPtmStr);
                        if (oldVarPtm != null) {
                            if (varPtm.priority > oldVarPtm.priority){
                                dstMap.put(varPtmStr, varPtm);
                            }
                        } else {
                            dstMap.put(varPtmStr, varPtm);
                        }

                        Set<Integer> allPoses = allMassAllPosesMap.get(mass);
                        if (allPoses != null) {
                            allPoses.add(absPos);
                        } else {
                            allPoses = new HashSet<>();
                            allPoses.add(absPos);
                            allMassAllPosesMap.put(varPtm.mass, allPoses);
                        }
                    }
                }

                if (!dstMap.isEmpty()) {
                    Map<Double, VarPtm> massVarPtmMap = new HashMap<>(dstMap.size());
                    for (VarPtm varPtm : dstMap.values()){
                        massVarPtmMap.put(varPtm.mass, varPtm);
                    }
                    pos_MassVarPtm_Map.put(absPos, massVarPtmMap);
                }
            }
        } // end for (int i = 0; i < partSeq.length(); ++i) {

        for (double mass : allMassAllPosesMap.keySet()){
            int times = allMassAllPosesMap.get(mass).size();
            if (times == 1) {
                absPos = Collections.min(allMassAllPosesMap.get(mass));

                Set<Double> massSet = oneTimeMassSetsMap.get(absPos);
                if (massSet != null) {
                    massSet.add(mass);
                } else {
                    massSet = new HashSet<>();
                    massSet.add(mass);
                    oneTimeMassSetsMap.put(absPos, massSet);
                }
            } else {
                Set<Integer> posComb = new HashSet<>(allMassAllPosesMap.get(mass));
                Set<Double> multiMassSet = posComb_multiMassSet_Map.get(posComb);
                if (multiMassSet != null) {
                    multiMassSet.add(mass);
                } else {
                    multiMassSet = new HashSet<>();
                    multiMassSet.add(mass);
                    for (int absPosTmp : posComb) {
                        if (oneTimeMassSetsMap.containsKey(absPosTmp)) {  // just think about PEPWD, pos 0 2 wont be in oneTimeMassSetsMap.keyset(), but P is in massWithMultiTime, so oneTimeMassSetsMap.get(pos) will fail.
                            multiMassSet.addAll(oneTimeMassSetsMap.get(absPosTmp));
                        }
                    }
                    posComb_multiMassSet_Map.put(posComb, multiMassSet);
                }
            }
        }
//        if (lszDebugScanNum.contains(scanNum)){
//            int a = 1;
//        }
    }

    public void prepareInfoNTerm(int scanNum, String partSeq,
                                 Map<Integer, Set<Double>> oneTimeMassSetsMap,
                                 Map<Set<Integer>, Set<Double>> posComb_multiMassSet_Map,
                                 Map<Integer, Map<Double, VarPtm>> pos_MassVarPtm_Map,
                                 Map<Double, Set<Integer>> allMassAllPosesMap,
                                 Map<Integer, Integer> yIdMinAbsPosMap,
                                 boolean couldBeProtN,
                                 int maxAbsNPos,
                                 int endRefPos,
                                 String protSeq) {

        int partSeqLen = partSeq.length();
        char aa;
        int absPos;
        for (int relPos = 0; relPos < partSeqLen; ++relPos) {
            aa = partSeq.charAt(relPos);
            absPos = relPos + endRefPos - partSeqLen;
            if (aaWithFixModSet.contains(aa) && absPos != 0) continue;  // if that pos has fix mod but is not N term, dont let it

            List<Byte> positionsToTry = new ArrayList<>(3);
            positionsToTry.add(ANYWHERE);
            if (yIdMinAbsPosMap.containsValue(absPos)) {
                positionsToTry.add(PEPN);
            }
            if (absPos == maxAbsNPos) {
                if (! nTermSpecific || absPos == 0 || protSeq.charAt(absPos-1) == 'K' || protSeq.charAt(absPos-1) == 'R') {
                    positionsToTry.add(PEPN);
                }
            }
            if (couldBeProtN && relPos == 0) {
                positionsToTry.add(PROTN);
            }

            if (aaAllVarPtmMap.containsKey(aa)) {
                Map<Byte, List<VarPtm>> allVarPtmMap = aaAllVarPtmMap.get(aa);
                Map<String, VarPtm> dstMap = new HashMap<>();
                for (Byte position : positionsToTry) {
                    if ( ! allVarPtmMap.containsKey(position)) continue;
                    for (VarPtm varPtm : allVarPtmMap.get(position)) {
                        double mass = varPtm.mass;
                        if (massTable.get(partSeq.charAt(relPos)) + mass < ms2Tolerance) continue;
                        String varPtmStr = varPtm.getStr();

                        VarPtm oldVarPtm = dstMap.get(varPtmStr);
                        if (oldVarPtm != null) {
                            if (varPtm.priority > oldVarPtm.priority){
                                dstMap.put(varPtmStr, varPtm);
                            }
                        } else {
                            dstMap.put(varPtmStr, varPtm);
                        }

                        Set<Integer> allPoses = allMassAllPosesMap.get(mass);
                        if (allPoses != null) {
                            allPoses.add(absPos);
                        } else {
                            allPoses = new HashSet<>();
                            allPoses.add(absPos);
                            allMassAllPosesMap.put(varPtm.mass, allPoses);
                        }
                    }
                }

                if (!dstMap.isEmpty()) {
                    Map<Double, VarPtm> massVarPtmMap = new HashMap<>();
                    for (VarPtm varPtm : dstMap.values()){
                        massVarPtmMap.put(varPtm.mass, varPtm);
                    }
                    pos_MassVarPtm_Map.put(absPos, massVarPtmMap);
                }
            }
        } // end for (int i = 0; i < partSeq.length(); ++i) {

        for (double mass : allMassAllPosesMap.keySet()){
            int times = allMassAllPosesMap.get(mass).size();
            if (times == 1) {
                absPos = Collections.min(allMassAllPosesMap.get(mass));

                Set<Double> massSet = oneTimeMassSetsMap.get(absPos);
                if (massSet != null) {
                    massSet.add(mass);
                } else {
                    massSet = new HashSet<>();
                    massSet.add(mass);
                    oneTimeMassSetsMap.put(absPos, massSet);
                }
            } else {
                Set<Integer> posComb = new HashSet<>(allMassAllPosesMap.get(mass));
                Set<Double> multiMassSet = posComb_multiMassSet_Map.get(posComb);
                if (multiMassSet != null) {
                    multiMassSet.add(mass);
                } else {
                    multiMassSet = new HashSet<>();
                    multiMassSet.add(mass);
                    for (int absPosTmp : posComb) {
                        if (oneTimeMassSetsMap.containsKey(absPosTmp)) {  // just think about PEPWD, pos 0 2 wont be in oneTimeMassSetsMap.keyset(), but P is in massWithMultiTime, so oneTimeMassSetsMap.get(pos) will fail.
                            multiMassSet.addAll(oneTimeMassSetsMap.get(absPosTmp));
                        }
                    }
                    posComb_multiMassSet_Map.put(posComb, multiMassSet);
                }
            }
        }
    }

    public void prepareInfoNTerm(int scanNum, String partSeq,
                                 Map<Integer, TreeMap<Double, VarPtm>> absPos_MassVarPtm_Map,
                                 Map<Integer, Integer> yIdMinAbsPosMap,
                                 boolean couldBeProtN,
                                 int optEndPosP1,
                                 int endRefPos,
                                 String protSeq,
                                 Map<Integer, List<Byte>> absPos_ptmPositions_Map) {

        int partSeqLen = partSeq.length();
        char aa;
        int absPos;
        for (int relPos = 0; relPos < partSeqLen; ++relPos) {
            aa = partSeq.charAt(relPos);
            absPos = relPos + endRefPos - partSeqLen;
            List<Byte> positionsToTry = new ArrayList<>(3);
            absPos_ptmPositions_Map.put(absPos, positionsToTry);
            if (aaWithFixModSet.contains(aa) && absPos != 0) continue;  // if that pos has fix mod but is not N term, dont let it

            positionsToTry.add(ANYWHERE);
            if (yIdMinAbsPosMap.containsValue(absPos)) {
                positionsToTry.add(PEPN);
            }
            if (absPos == optEndPosP1) {
                positionsToTry.add(PEPN);
            }
            if (couldBeProtN && relPos == 0) {
                positionsToTry.add(PROTN);
            }

            if (aaAllVarPtmMap.containsKey(aa)) {
                Map<Byte, List<VarPtm>> allVarPtmMap = aaAllVarPtmMap.get(aa);
                Map<String, VarPtm> dstMap = new HashMap<>();
                for (Byte position : positionsToTry) {
                    if ( ! allVarPtmMap.containsKey(position)) continue;
                    for (VarPtm varPtm : allVarPtmMap.get(position)) {
                        double mass = varPtm.mass;
                        if (massTable.get(partSeq.charAt(relPos)) + mass < ms2Tolerance) continue;
                        String varPtmStr = varPtm.getStr();

                        VarPtm oldVarPtm = dstMap.get(varPtmStr);
                        if (oldVarPtm != null) {
                            if (varPtm.priority > oldVarPtm.priority){
                                dstMap.put(varPtmStr, varPtm);
                            }
                        } else {
                            dstMap.put(varPtmStr, varPtm);
                        }

                    }
                }

                if (!dstMap.isEmpty()) {
                    TreeMap<Double, VarPtm> massVarPtmMap = new TreeMap<>();
                    for (VarPtm varPtm : dstMap.values()){
                        massVarPtmMap.put(varPtm.mass, varPtm);
                    }
                    absPos_MassVarPtm_Map.put(absPos, massVarPtmMap);
                }
            }
        } // end for (int i = 0; i < partSeq.length(); ++i) {

    }
}
