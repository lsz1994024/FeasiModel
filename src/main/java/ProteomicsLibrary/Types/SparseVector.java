package ProteomicsLibrary.Types;

import java.util.*;

public class SparseVector {

    private Map<Integer, Double> sparseVector = new HashMap<>();


    public SparseVector() {}

    public void add(int i, double v) {
        if (Math.abs(v) > 1e-6) {
            if (sparseVector.containsKey(i)) {
                sparseVector.put(i, sparseVector.get(i) + v);
            } else {
                sparseVector.put(i, v);
            }
        }
    }

    public void put(int i, double v) {
        if (Math.abs(v) > 1e-6) {
            sparseVector.put(i, v);
        }
    }

    public double get(int i) {
        if (sparseVector.containsKey(i)) {
            return sparseVector.get(i);
        } else {
            return 0;
        }
    }
    public Double[] getValues() {
        return sparseVector.values().toArray(new Double[0]);
    }
    public boolean isEmpty() {
        return sparseVector.isEmpty();
    }

}
