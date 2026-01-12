import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.ArrayList;
public class SchFactAna {
    private static SchFactAna fact = new SchFactAna();
    public static SchFactAna getInstance() { return fact; }
    public double[][] primFactMethod(double[][] xij) {
        FactAna_0 fact_0 = new PrimFactMethod(xij);

        return fact_0.factLoading();
    }
    public double[][] maxLikeMethod(double[][] xij) {
        FactAna_0 fact_0 = new MaxLikeMethod(xij);

        return fact_0.factLoading();
    }
    public Contribution2[] contribution(double[][] factlds) {
        FactAna_1 fact_1 = new FactAna_1(factlds);

        return fact_1.contribution();
    }
    /*********************************/
    /* interface define              */
    /*********************************/
    /*********************************/
    /* class define                  */
    /*********************************/
    // 寄与率・累積寄与率
    public class Contribution2 {
        private double cr = 0.0;
        private double ccr = 0.0;
        public Contribution2(double cr, double ccr) {
            this.cr = cr;
            this.ccr = ccr;
        }
        public double getCr()      { return this.cr; }
        public double getCcr()     { return this.ccr; }
    }

    private abstract class FactAna_0 {
        private RealMatrix matrixAT = null;      // 元データ
        private int rowN = 0;
        private int colN = 0;
        abstract double[][] factExtract(double[][] corr);
        public FactAna_0(RealMatrix matrixAT) {
            this.matrixAT = matrixAT;
            this.rowN = matrixAT.getRow(0).length;
            this.colN = matrixAT.getColumn(0).length;
        }
        public double[][] factLoading() {
            double[][] corr = calcCorrMatrix(matrixAT).getData();  // 相関行列
            double[][] retFact = factExtract(corr);

            if (1 < retFact[0].length) {
                retFact = factRol(retFact);
            }
            return retFact;
        }
        // 因子回転
        private double[][] factRol(double[][] retFact) {
            RealMatrix matrixB = MatrixUtils.createRealMatrix(retFact);
            RealMatrix matrixR = calcCorrMatrix(matrixB.transpose());

            return matrixB.multiply(matrixR).getData();
        }

        // 相関行列作成
        private RealMatrix calcCorrMatrix(RealMatrix matrixWT) {
            int colR = matrixWT.getColumn(0).length;
            PearsonsCorrelation corel = new PearsonsCorrelation();
            double[][] matrix = new double[colR][colR];
            
            for(int i = 0; i < colR; i++) {
                for(int j = 0; j < colR; j++) {
                    double[] xArray = matrixWT.getRow(i);
                    double[] yArray = matrixWT.getRow(j);

                    matrix[i][j] = corel.correlation(xArray, yArray);
                }
            }
            return MatrixUtils.createRealMatrix(matrix);
        }
    }

    private class FactAna_1 {
        private RealMatrix matrixAT = null;      // 元データ
        private int rowN = 0;
        private int colN = 0;
        private double sum = 0.0;
        public FactAna_1(double[][] xij) {
            this.matrixAT = MatrixUtils.createRealMatrix(xij).transpose();
            this.rowN = matrixAT.getRow(0).length;
            this.colN = matrixAT.getColumn(0).length;
        }
        // 寄与率・累積寄与率
        public Contribution2[] contribution() {
            double[] lambas = calcLamba();
            Contribution2[] contris = new Contribution2[colN];
            double ccr = 0.0;

            for(int i = 0; i < colN; i++) {
                double cr = lambas[i] / sum;

                ccr += cr;
                contris[i] = new Contribution2(cr, ccr);
            }
            return contris;
        }
        private double[] calcLamba() {
            double[] lamba = new double[colN];

            for(int i = 0; i < colN; i++) {
                double[] factld = matrixAT.getRow(i);
                double s = 0.0;

                for(int j = 0; j < factld.length; j++) {
                    s += factld[j] * factld[j];
                }
                lamba[i] = s / factld.length;
                sum += lamba[i];
            }
            return lamba;
        }
    }
    // 主因子法
    private class PrimFactMethod extends FactAna_0 {
        int colN = 0;
        public PrimFactMethod(double[][] xij) {
            super(MatrixUtils.createRealMatrix(xij).transpose());
        }
        // 因子抽出
        double[][] factExtract(double[][] corr) {
            colN = corr.length;
            double[][] retFact = null;
            double[][] wkFact = new double[colN][colN];
            boolean factZ = true;
            int factN = 0;             // 因子数

            do {
                factZ = true;
                Eigen[] eds = eigen(corr);
                factN = eds.length;
                int j = 0;
                
                for(Eigen ed : eds) {
                    double lambSq = Math.sqrt(ed.getEdVal());
                    double[] edVec = ed.getEdVec();

                    for(int i = 0; i < colN; i++) {
                        wkFact[i][j] = lambSq * edVec[i];
                    }
                    j++;
                }
                // 収束判定
                for(int i = 0; i < colN; i++) {
                    double wk = 0.0;

                    for(int p = 0; p < factN; p++) {
                        wk += wkFact[i][p] * wkFact[i][p];
                    }
                    if (Math.abs(wk - corr[i][i]) > 0.0005) {
                        factZ = false;
                    }
                    corr[i][i] = wk;   // 相関行列更新
                }
            } while(factZ != true);
            retFact = new double[colN][factN];
            for(int i = 0; i < colN; i++) {
                retFact[i] = Arrays.copyOf(wkFact[i], factN);
            }
            return retFact;
        }
        // 小数点以下四捨五入
        private double roundEx(double val, int expo) {
            double result = Math.pow(10, expo);

            return Math.round(val*result) / result;
        }
        // 小数点以下切り捨て
        private double floorEx(double val, int expo) {
            double result = Math.pow(10, expo);

            return Math.floor(val*result) / result;
        }
        // 固有値・固有ベクトル
        private Eigen[] eigen(double[][] corr) {
            List<Eigen> arrEds = new ArrayList<Eigen>();
            // 固有値計算
            RealMatrix matrixB = MatrixUtils.createRealMatrix(corr);
            EigenDecomposition ed = new EigenDecomposition(matrixB);
            double[] eigens = ed.getRealEigenvalues();
            // 固有値を求め、それの固有ベクトルを求める
            for(int i = 0; i < eigens.length; i++) {
                double[] evs = ed.getEigenvector(i).toArray();

                if (1 <= eigens[i]) {
                    arrEds.add(new Eigen(eigens[i], evs));
                }
            }
            arrEds.sort(Comparator.comparing(Eigen::getEdVal).reversed());
            return arrEds.toArray(new Eigen[arrEds.size()]);
        }
    }
    // 最尤法
    private class MaxLikeMethod extends FactAna_0 {
        private double[][] xij = null;
        public MaxLikeMethod(double[][] xij) {
            super(MatrixUtils.createRealMatrix(xij));
            this.xij = xij;
        }
        // 因子抽出
        double[][] factExtract(double[][] corr) {
            Eigen[] eds = eigen(corr);
            return corr;
        }
        // 固有値・固有ベクトル
        private Eigen[] eigen(double[][] corr) {
            List<Eigen> arrEds = new ArrayList<Eigen>();
            // 固有値計算
            RealMatrix matrixB = MatrixUtils.createRealMatrix(corr);
            EigenDecomposition ed = new EigenDecomposition(matrixB);
            double[] eigens = ed.getRealEigenvalues();
            // 固有値を求め、それの固有ベクトルを求める
            for(int i = 0; i < eigens.length; i++) {
                double[] evs = ed.getEigenvector(i).toArray();

                    arrEds.add(new Eigen(eigens[i], evs));
            }
            return arrEds.toArray(new Eigen[arrEds.size()]);
        }
    }
}

