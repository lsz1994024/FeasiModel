package proteomics.Segment;


import ProteomicsLibrary.DbTool;
import ProteomicsLibrary.MassTool;
import ProteomicsLibrary.Types.*;
import org.apache.commons.math3.util.Pair;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Types.*;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static proteomics.PIPI.MIN_PEAK_SUM_INFER_AA;

public class InferSegment {
    public static final byte N_TAG = -1;
    public static final byte NON_NC_TAG = 0;
    public static final byte C_TAG = 1;
    private final double ms2Tolerance;
    private TreeMap<Segment, Integer> aaVectorTemplate = new TreeMap<>();
    private Map<Double, String> augedMassAaMap = new HashMap<>(35, 1);// augmented amino acid map. Normal aa plus aa with mod
    private MassTool massTool;
    public Map<String, Double> PTM_mass_Map = new HashMap<>(); // used for generating a correct tag
    public InferSegment(MassTool massTool, Map<String, String> parameterMap, Map<Character, Double> fixModMap) throws Exception {
        PTM_mass_Map.put("Phospho",79.966);//4
        PTM_mass_Map.put("Oxidation",15.995); //12
        PTM_mass_Map.put("Carbamidomethyl",57.021);
        PTM_mass_Map.put("Fluoro",17.991);//10
        PTM_mass_Map.put("Carboxyethyl", 72.021);//7
        PTM_mass_Map.put("Chlorination",33.961 );//13
        PTM_mass_Map.put("Malonyl",86.000);//5
        PTM_mass_Map.put("Methylamine",13.032);//14
        PTM_mass_Map.put("Nitro", 44.985);//11
        PTM_mass_Map.put("Sulfide",31.972);//8
        PTM_mass_Map.put("Deamidated",0.984);//15
        PTM_mass_Map.put("Nitrosyl",28.990 );//1
        PTM_mass_Map.put("Decarboxylation",-30.011);//6
        PTM_mass_Map.put("Deoxyhypusine", 71.073);//2
        PTM_mass_Map.put("Dioxidation",31.990);//17
        PTM_mass_Map.put("Ethanolyl", 44.026);//9
        PTM_mass_Map.put("Carbonyl",13.979);//3


        this.massTool = massTool;
        this.ms2Tolerance = Double.valueOf(parameterMap.get("ms2_tolerance"));
        Map<Character, Double> massTable = massTool.getMassTable();

        char[] standardAaArray = new char[]{'G', 'A', 'S', 'P', 'V', 'T', 'C', 'L', 'N', 'D', 'Q', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W'};

        Map<Double, Character> oriMassAaMap = new HashMap<>(25, 1);
        for (char aa : standardAaArray) {
            // # = I/L.
            oriMassAaMap.put(massTable.get(aa), aa);
        }

        Character[] aaArray = oriMassAaMap.values().toArray(new Character[0]);

        for (char aa1 : aaArray) {
            for (char aa2 : aaArray) {
                for (char aa3 : aaArray) {
                    aaVectorTemplate.put(new Segment(String.format(Locale.US, "%c%c%c", aa1, aa2, aa3)), 0);
                }
            }
        }

        int idx = 0;
        for (Segment segment : aaVectorTemplate.keySet()) {
            aaVectorTemplate.put(segment, idx);
            ++idx;
        }

        // generate a mass aa map containing modified amino acid
        for (double k : oriMassAaMap.keySet()) {
            augedMassAaMap.put(k, oriMassAaMap.get(k).toString());
        }
    }

    public SparseBooleanVector generateSegmentBooleanVector(String peptide) {

        String normalizedPeptide = normalizeSequence(DbTool.getSequenceOnly(peptide));

        Set<Integer> tempSet = new HashSet<>(DbTool.getSequenceOnly(peptide).length() + 1, 1);
        for (int i = 0; i <= normalizedPeptide.length() - 3; ++i) {
            tempSet.add(aaVectorTemplate.get(new Segment(normalizedPeptide.substring(i, i + 3))));
        }
        return new SparseBooleanVector(tempSet);
    }
    public static String normalizeSequence(String seq) {
        return seq.replaceAll("I", "L");
    }

    public ExpTag getOneCorrectTag( List<String> truthStrList)  {
        String freeSeq = truthStrList.get(0).replace('I','L');
        String[] mods = truthStrList.get(2).split(";");
        String[] mod_sites = truthStrList.get(3).split(";");

        int modNum = mods.length;
        int seqLen = freeSeq.length();

        int minFreePos = 0;
        int tagLPos;
        int tagLen;
        Peptide peptide = new Peptide(freeSeq, false, massTool);

        int maxFreeGapLen = 0;
        if (! mods[0].isEmpty()){// mods[0].isEmpty() means mods = [""], i.e. no mods for this scan
            int lastPos = -1;
            PosMassMap posMassMap = new PosMassMap();
            for (int ptmId = 0; ptmId < modNum; ptmId++) {
                double mass = PTM_mass_Map.get(mods[ptmId].split("@")[0]);
                int pos = Integer.valueOf(mod_sites[ptmId]) - 1;
                posMassMap.put(pos, mass);
                if (pos-lastPos-1 > maxFreeGapLen) {
                    maxFreeGapLen = pos-lastPos - 1;
                    minFreePos = lastPos+1;
                }
                lastPos = pos;
            }
            if (seqLen-lastPos-1 > maxFreeGapLen) { // last gap
                maxFreeGapLen = seqLen-lastPos - 1;
                minFreePos = lastPos+1;
            }
            peptide.setVarPTM(posMassMap);
        } else {
            maxFreeGapLen = seqLen - minFreePos - 1;
        }
        double[][] ionMatrix = peptide.getIonMatrixNow();
        tagLen = (int) (maxFreeGapLen/1.2); // must be an int
//        tagLen = maxFreeGapLen - 3; // must be an int
        if (maxFreeGapLen <= 2) {
            return null;
        } else if (maxFreeGapLen <= 5) {
            tagLen = 3;
        }
        tagLPos = minFreePos + (maxFreeGapLen-tagLen)/2; // must be an int

        if (minFreePos == 0) { // possible to be Ntag
            boolean chance = ((int) (ionMatrix[0][ionMatrix[0].length-1]) % 2) == 0;
            if (chance) {
                tagLPos = 0;
            }
        }
        if (minFreePos + maxFreeGapLen == seqLen) { //possible to be Ctag
            boolean chance = ((int) (ionMatrix[1][2]) % 2) == 0;
            if (chance) {
                tagLPos = seqLen - tagLen;
            }
        }

        List<ExpAa> expAaList = new ArrayList<>();
        for (int aaId = tagLPos; aaId < tagLPos+tagLen; aaId++) {
            int isNorC = 0;
            if (aaId == 0) {
                isNorC = -1;
            } else if (aaId == seqLen-1) {
                isNorC = 1;
            }
            double headMz;
            if (aaId == 0) {
                headMz = MassTool.PROTON;
            } else {
                headMz = ionMatrix[0][aaId-1];
            }
            double tailMz = ionMatrix[0][aaId];
            ExpAa aa = new ExpAa(freeSeq.substring(aaId, aaId+1), freeSeq.charAt(aaId), headMz, tailMz, 1, 1, isNorC);
            expAaList.add(aa);
        }
        if (expAaList.isEmpty()) return null;
        ExpTag resTag = new ExpTag(expAaList);
        if (expAaList.get(0).isNorC == -1) {
            resTag.isNorC = -1;
        } else if (expAaList.get(expAaList.size()-1).isNorC == 1) {
            resTag.isNorC = 1;
        }
        return resTag;
    }


    public TreeMap<Double, Double> addVirtualPeaks(double precursorMass, TreeMap<Double, Double> plMap) {
        double totalMass = precursorMass + 2 * MassTool.PROTON;
        TreeMap<Double, Double> finalPlMap = new TreeMap<>();
        for (double mz : plMap.keySet()) {
            finalPlMap.put(mz, plMap.get(mz));
        }
        for (double mz : plMap.keySet()) {
            double anotherMz = totalMass - mz;
            double leftMz = anotherMz - ms2Tolerance;
            double rightMz = anotherMz + ms2Tolerance;
            NavigableMap<Double, Double> temp = null;
            try {
                temp = plMap.subMap(leftMz, true, rightMz, true);
            } catch (IllegalArgumentException ex) {}

        }

        // Add two virtual peak. Because we have convert all y-ions to b-ions.
        finalPlMap.put(MassTool.PROTON, 1d);
        double cTermMz = precursorMass - massTool.H2O + MassTool.PROTON;
        double leftMz = cTermMz - ms2Tolerance;
        double rightMz = cTermMz + ms2Tolerance;
        NavigableMap<Double, Double> temp = null;
        try {
            temp = plMap.subMap(leftMz, true, rightMz, true);
        } catch (IllegalArgumentException ex) {}
        if ((temp == null) || (temp.isEmpty())) {
            finalPlMap.put(cTermMz, 1d);
        }
        finalPlMap.put(MassTool.PROTON + massTool.H2O, 1d);

        return finalPlMap;
    }
}
