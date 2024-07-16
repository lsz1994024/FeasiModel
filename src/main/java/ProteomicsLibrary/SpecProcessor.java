package ProteomicsLibrary;

import ProteomicsLibrary.Types.SparseVector;

import java.util.*;

public class SpecProcessor {

    private static final double defaultIntensity = 1; // DO NOT change. Otherwise, change the whole project accordingly.
    private static final double removePrecursorPeakTolerance = 1.5; // this equals the isolation window.
    private final MassTool massTool;

    public SpecProcessor(MassTool massTool) {
        this.massTool = massTool;
    }

    public TreeMap<Double, Double> preSpectrumTopNStyleWithChargeLimit (Map<Double, Double> inputPL, double precursorMass, int precursorCharge, double minClear, double maxClear, int topN, double ms2Tolerance) {

        TreeMap<Double, Double> outputPL = removeCertainPeaks(inputPL, precursorMass, precursorCharge, minClear, maxClear);
        List<Map.Entry<Double, Double>> plList = new ArrayList<>(outputPL.entrySet());
        Collections.sort(plList, Comparator.comparing(o -> o.getValue(), Comparator.reverseOrder()));
        plList = plList.subList(0, Math.min(plList.size(), 600));
        Collections.sort(plList, Comparator.comparing(o -> o.getKey()));

        TreeMap<Double, Double> newPlMap = new TreeMap<>();

        for (Map.Entry<Double, Double> peak : plList ) {
            newPlMap.put(peak.getKey(), peak.getValue());
        }
        if (newPlMap.subMap(0d, precursorMass).isEmpty()) {
            return new TreeMap<>();
        } else {
            return topNStyleNormalization(sqrtPL(new TreeMap<>(newPlMap.subMap(0d, precursorMass))), topN); //todo this deletes the precursorMass peak, ie the b-end ion and y-end ion
        }
    }

    public SparseVector digitizePL(TreeMap<Double, Double> plMap) {
        SparseVector digitizedPL = new SparseVector();
        for (double mz : plMap.keySet()) {
            int idx = massTool.mzToBin(mz);
            if (Math.abs(plMap.get(mz)) > 1e-6) {
                digitizedPL.put(idx, Math.max(plMap.get(mz), digitizedPL.get(idx)));
            }
        }
        return digitizedPL;
    }

    public static TreeMap<Double, Double> topNStyleNormalization(TreeMap<Double, Double> inputPL, int localTopN) {
        if (inputPL.isEmpty()) {
            return new TreeMap<>();
        } else {
            // select top N in each 100 Da
            TreeMap<Double, Double> outputPL = new TreeMap<>();
            double minMz = inputPL.firstKey();
            double maxMz = inputPL.lastKey();
            double leftMz = minMz;
            while (leftMz < maxMz) {
                // find the max intensity in each window
                double rightMz = Math.min(leftMz + 100, maxMz);
                NavigableMap<Double, Double> subMap;
                if (rightMz < maxMz) {
                    subMap = inputPL.subMap(leftMz, true, rightMz, false);
                } else {
                    subMap = inputPL.subMap(leftMz, true, rightMz, true);
                }
                if (!subMap.isEmpty()) {
                    Double[] intensityArray = subMap.values().toArray(new Double[0]);
                    Arrays.sort(intensityArray, Comparator.reverseOrder());
                    double temp1 = defaultIntensity / intensityArray[0];
                    double temp2 = subMap.size() > localTopN ? intensityArray[localTopN] : 0;
                    for (double mz : subMap.keySet()) {
                        if (subMap.get(mz) > temp2) {
                            outputPL.put(mz, subMap.get(mz) * temp1);
                        }
                    }
                }
                leftMz = rightMz;
            }

            return outputPL;
        }
    }

    private static TreeMap<Double, Double> removeCertainPeaks(Map<Double, Double> peakMap, double precursorMass, int precursorCharge, double minClear, double maxClear) {
        TreeMap<Double, Double> mzIntensityMap = new TreeMap<>();
        double precursorMz = precursorMass / precursorCharge + MassTool.PROTON;
        for (double mz : peakMap.keySet()) {
            if (((mz < minClear) || (mz > maxClear)) && (mz > 50)) {
                if ((peakMap.get(mz) > 1e-6) && (Math.abs(peakMap.get(mz) - precursorMz) > removePrecursorPeakTolerance)) {
                    mzIntensityMap.put(mz, peakMap.get(mz));
                }
            }
        }

        return mzIntensityMap;
    }
    private static TreeMap<Double, Double> sqrtPL(TreeMap<Double, Double> plMap) {
        // sqrt the intensity and find the highest intensity.
        TreeMap<Double, Double> sqrtPlMap = new TreeMap<>();
        for (double mz : plMap.keySet()) {
            if (plMap.get(mz) > 1e-6) {
                sqrtPlMap.put(mz, Math.sqrt(plMap.get(mz)));
            }
        }
        return sqrtPlMap;
    }
}
