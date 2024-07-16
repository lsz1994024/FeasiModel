package proteomics.Types;
import java.util.List;

public class ExpTag implements Comparable<ExpTag> {
    public int isNorC = 0; //N -1, none 0, C 1;
    private final ExpAa[] expAaArray;
    private int hashCode;
    private final double totalIntensity;
    private final String freeAaString;
    private final String ptmAaString;
    public List<ExpAa> expAaList = null;
    public double[] intensityArray;
    public ExpTag(ExpAa aa1, ExpAa aa2, ExpAa aa3) {
        expAaArray = new ExpAa[]{aa1, aa2, aa3};
        String toString = expAaArray[0].toString() + "-" + expAaArray[1].toString() + "-" + expAaArray[2].toString();
        hashCode = toString.hashCode();

        StringBuilder sbFreeSeq = new StringBuilder(expAaArray.length);
        StringBuilder sbPtmSeq = new StringBuilder(2* expAaArray.length - 1);
        for (ExpAa aa : expAaArray) {
            sbFreeSeq.append(aa.getPtmFreeAA());
            sbPtmSeq.append(aa.getAA());
        }
        freeAaString = sbFreeSeq.toString();
        ptmAaString = sbPtmSeq.toString();

        double intensity = expAaArray[0].getHeadIntensity();
        for (ExpAa aa : expAaArray) {
            intensity += aa.getTailIntensity();
        }
        totalIntensity = intensity;
    }
    public ExpTag(ExpAa aa1, ExpAa aa2, ExpAa aa3, ExpAa aa4) {
        expAaArray = new ExpAa[]{aa1, aa2, aa3, aa4};
        String toString = expAaArray[0].toString() + "-" + expAaArray[1].toString() + "-" + expAaArray[2].toString()+ "-" + expAaArray[3].toString();
        hashCode = toString.hashCode();

        StringBuilder sbFreeSeq = new StringBuilder(expAaArray.length);
        StringBuilder sbPtmSeq = new StringBuilder(2* expAaArray.length - 1);
        for (ExpAa aa : expAaArray) {
            sbFreeSeq.append(aa.getPtmFreeAA());
            sbPtmSeq.append(aa.getAA());
        }
        freeAaString = sbFreeSeq.toString();
        ptmAaString = sbPtmSeq.toString();

        double intensity = expAaArray[0].getHeadIntensity();
        for (ExpAa aa : expAaArray) {
            intensity += aa.getTailIntensity();
        }
        totalIntensity = intensity;
    }
    public ExpTag(List<ExpAa> expAaList) {
        isNorC = expAaList.get(0).isNorC;
        intensityArray = new double[expAaList.size()+1];
        this.expAaList = expAaList;
        intensityArray[0] = expAaList.get(0).getHeadIntensity();
        for (int i = 0; i < expAaList.size(); i++) {
            intensityArray[i+1] = expAaList.get(i).getTailIntensity();
        }
        expAaArray = expAaList.toArray(new ExpAa[0]);
        String toString = expAaArray[0].toString();
        for (int i = 1; i < expAaArray.length; i++){
            toString += "-" + expAaArray[i].toString();
        }
        hashCode = toString.hashCode();
        StringBuilder sbFreeSeq = new StringBuilder(expAaArray.length);
        StringBuilder sbPtmSeq = new StringBuilder(2* expAaArray.length - 1);
        for (ExpAa aa : expAaArray) {
            sbFreeSeq.append(aa.getPtmFreeAA());
            sbPtmSeq.append(aa.getAA());
        }
        freeAaString = sbFreeSeq.toString();
        ptmAaString = sbPtmSeq.toString();

        double intensity = expAaArray[0].getHeadIntensity();
        for (ExpAa aa : expAaArray) {
            intensity += aa.getTailIntensity();
        }
        totalIntensity = intensity;
    }

    public int hashCode() {
        return hashCode;
    }

    public boolean equals(Object other) {
        return (other instanceof ExpTag) && (this.hashCode() == other.hashCode());
    }
    public int compareTo(ExpTag other) {
        return Double.compare(expAaArray[0].getHeadLocation(), other.expAaArray[0].getHeadLocation());
    }

    public String getFreeAaString() {
        return freeAaString;
    }
    public String getPtmAaString() {
        return ptmAaString;
    }

    public String toString() {
        return isNorC+","+ ptmAaString ;
    }

    public double getTotalIntensity() {
        return totalIntensity;
    }

    public double getHeadLocation() {
        return expAaArray[0].getHeadLocation();
    }

    public double getTailLocation() {
        return expAaArray[expAaArray.length - 1].getTailLocation();
    }

    public ExpTag clone() throws CloneNotSupportedException {
        super.clone();
        if (expAaArray.length == 3 ){
            return new ExpTag(expAaArray[0].clone(), expAaArray[1].clone(), expAaArray[2].clone());
        }else {
            return new ExpTag(expAaArray[0].clone(), expAaArray[1].clone(), expAaArray[2].clone(), expAaArray[3].clone());
        }
    }

    public int size() {
        return expAaArray.length;
    }

    public ExpAa get(int i) {
        return expAaArray[i];
    }

}
