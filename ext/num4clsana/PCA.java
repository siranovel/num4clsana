import org.apache.commons.math3.stat.correlation.Covariance;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import org.apache.commons.math3.stat.StatUtils;

import java.util.Arrays;
import java.util.Comparator;
public class PCA {
    private static PCA pca = new PCA();
    public static PCA getInstance() { return pca; }
    public Eigen[] eigen(double[][] xij) {
        PCA_1 pca_1 = new PCA_1(xij);

        return pca_1.eigen();
    }
    public Contribution[] contribution(Eigen[] eds, double[][] xij) {
        PCA_1 pca_1 = new PCA_1(xij);

        return pca_1.contribution(eds);
    }
    public Score[] score(Eigen[] eds, double[][] xij) {
        PCA_1 pca_1 = new PCA_1(xij);

        return pca_1.score(eds);
    }
    /*********************************/
    /* interface define              */
    /*********************************/
    /*********************************/
    /* class define                  */
    /*********************************/
    // 寄与率・累積寄与率
    public class Contribution {
        private double edVal = 0.0;
        private double cr = 0.0;
        private double ccr = 0.0;
        public Contribution(double edVal, double cr, double ccr) {
            this.edVal = edVal;
            this.cr = cr;
            this.ccr = ccr;
        }
        public double getEdVal()   { return this.edVal; }
        public double getCr()      { return this.cr; }
        public double getCcr()     { return this.ccr; }
    }
    public class Score {
        private double edVal = 0.0;
        private double[] score = null;
        public Score(double edVal, double[] score) {
            this.edVal = edVal;
            this.score = score;
        }
        public double getEdVal() { return this.edVal; }
        public double[] getScore() { return this.score; }
    }
    // 主成分分析
    private class PCA_1 {
        private RealMatrix matrixAT = null;
        private int rowN = 0;
        private int colN = 0;
        public  PCA_1(double[][] xij) {
            this.matrixAT = MatrixUtils.createRealMatrix(xij).transpose();
            this.rowN = matrixAT.getRow(0).length;
            this.colN = matrixAT.getColumn(0).length;
        }
        // 固有値・固有ベクトル
        public Eigen[] eigen() {
            // 分散行列作成
            RealMatrix matrixA = calcCovatrianceMatrix();
            // 固有値・固有ベクトル計算
            Eigen[] eds = calcEigenVector(matrixA);
            // 固有値を小さい順にソート
            Arrays.sort(eds, Comparator.comparing(Eigen::getEdVal));
            return eds;
        }
        private RealMatrix calcCovatrianceMatrix(){
            Covariance corel = new Covariance();
            double[][] omg = new double[colN][colN];

            for(int i = 0; i < colN; i++) {
                for(int j = 0; j < colN; j++) {
                    double[] xArray = matrixAT.getRow(i);
                    double[] yArray = matrixAT.getRow(j);

                    omg[i][j] = corel.covariance(xArray, yArray);
                }
            }
            return MatrixUtils.createRealMatrix(omg);
        }
        private Eigen[] calcEigenVector(RealMatrix matrixA) {
            // 固有値計算
            EigenDecomposition ed = new EigenDecomposition(matrixA);
            double[] eigens = ed.getRealEigenvalues();
            Eigen[] eds = new Eigen[eigens.length];
            // 固有値を求め、それの固有ベクトルを求める
            for(int i = 0; i < eigens.length; i++) {
                double[] evs = ed.getEigenvector(i).toArray();
                eds[i] = new Eigen(eigens[i], evs);
            }
            return eds;
        }
        // 累積寄与率・寄与率
        public Contribution[] contribution(Eigen[] eds) {
            Contribution[] contris = new Contribution[eds.length];
            double info = calcBeforeInfo();
            int i = 0;
            double ccr = 0.0;

            for(Eigen ed : eds) {
                double edVal = ed.getEdVal();
                double cr = (info - edVal) / info;                

                ccr += cr;
                contris[i] = new Contribution(edVal, cr, ccr);
                i++;
            }
            return contris;
        }
        private double calcBeforeInfo() {
            double info = 0.0;
            
            for(int i = 0; i < colN; i++) {
                double[] dt = matrixAT.getRow(i);

                info += StatUtils.variance(dt);
            }
            return info;
        }
        // 主成分得点
        public Score[] score(Eigen[] eds) {
            Score[] scs = new Score[eds.length];
            double[] means = calcMean();
            int i = 0;

            for(Eigen ed : eds) {
                double edVal = ed.getEdVal();
                double[] edVec = ed.getEdVec();

                scs[i] = new Score(edVal, calcScore(edVec, means));
                i++;
            }
            return scs;
        }
        // 平均を求める
        private double[] calcMean() {
            double[] means = new double[colN];

            for(int i= 0; i < means.length; i++) {
                double[] dt = matrixAT.getRow(i);

                means[i] = StatUtils.mean(dt);
            }
            return means;
        }
        // 主成分得点を計算する
        private double[] calcScore(double[] edVec, double[] means) {
            double[] scores = new double[rowN];

            for(int i = 0; i < edVec.length; i++) {
                double[] dt = matrixAT.getRow(i);

                for(int j = 0; j < dt.length; j++) {
                    scores[j] += edVec[i] * (dt[j] - means[i]);
                }
            }
            return scores;
        }
    }
}

