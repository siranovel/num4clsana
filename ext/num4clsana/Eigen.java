public class Eigen {
    private double edVal= 0.0;     // 固有値
    private double[] edVec = null; // 固有ベクトル
    public Eigen(double edVal, double[] edVec) {
        this.edVal = edVal;
        this.edVec = edVec;
    }
    public double getEdVal()   { return this.edVal; }
    public double[] getEdVec() { return this.edVec; }
}

