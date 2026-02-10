import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.Covariance;
import org.apache.commons.math3.stat.StatUtils;

public class DiscAna {
    private static DiscAna da = new DiscAna();
    public static DiscAna getInstance() { return da;}
    public double[] line_disc(double[][] xa, double[][] xb){
        DiscAna_0 da0 = new DiscAna_0();

        return da0.line_disc(xa, xb);
    }
    public double[] score(double[] an, double[][] gx){
        DiscAna_0 da0 = new DiscAna_0();

        return da0.score(an, gx);
    }
    /*********************************/
    /* interface define              */
    /*********************************/
    /*********************************/
    /* class define                  */
    /*********************************/
    private class DiscAna_0 {
        private int n1 = 0;
        private int n2 = 0;
        // 線形型判別
        public double[] line_disc(double[][] xa, double[][] xb){
            // 各グループの個数
            RealMatrix matrixXa = MatrixUtils.createRealMatrix(xa);
            n1 = matrixXa.getColumn(0).length;
            RealMatrix matrixXb = MatrixUtils.createRealMatrix(xb);
            n2 = matrixXb.getColumn(0).length;

            // 各グループの平均を計算
            double[] meanG1 = calcMean(matrixXa);
            double[] meanG2 = calcMean(matrixXb);
            double[] meanG = new double[] {
                                meanG1[0] - meanG2[0],
                                meanG1[1] - meanG2[1],
                             };
            // 各グループの分散共分散行列を計算
            double[][] corrG1 = calcCorrMatrix(matrixXa);
            double[][] corrG2 = calcCorrMatrix(matrixXb);
            RealMatrix corrG = addMatrix(corrG1, corrG2);

            // 連立方程式を解く
            double[] wb = simEqu(corrG, meanG);
            double[] bn = new double[3];

            bn[0] = calcA0(wb, meanG1, meanG2) / wb[1];
            bn[1] = wb[0] / wb[1];
            bn[2] = wb[1] / wb[1];

            return bn;
        }
        private double[] calcMean(RealMatrix matrixX) {
            double[] m = new double[2];

            for(int i = 0; i < 2; i++) {
                double[] dt = matrixX.getColumn(i);

                m[i] = StatUtils.mean(dt);
            }
            return m;
        }
        private double[][] calcCorrMatrix(RealMatrix matrixX) {
            Covariance corel = new Covariance();
            double[][] corr = new double[2][2];

            for(int i =0; i < 2; i++) {
                for(int j = 0; j < 2; j++) {
                    double[] xArray = matrixX.getColumn(i);
                    double[] yArray = matrixX.getColumn(j);
                    
                    corr[i][j] = corel.covariance(xArray, yArray);
                }
            }
            return corr;
        }
        private RealMatrix addMatrix(double[][] corrG1,double[][] corrG2) {
            double[][] matrixG = new double[2][2];

            for(int i = 0; i < 2; i++) {
                for(int j = 0; j < 2; j++) {
                    matrixG[i][j] = 
                        ((n1 - 1) * corrG1[i][j] + (n2 - 1) * corrG2[i][j]) 
                      / (n1 + n2 - 2);
                }
            }
            return MatrixUtils.createRealMatrix(matrixG);
        }
        private double[] simEqu(RealMatrix corrG, double[] meanG) {
            RealMatrix corrGi = MatrixUtils.blockInverse(corrG, 0);
            double[] matrixRet = corrGi.operate(meanG);

            return matrixRet;
        }
        private double calcA0(double[] wb, double[] meanG1,double[] meanG2) {
            double ret = 0.0;

            for(int i = 0; i < 2; i++) {
                ret += wb[i] * (n1 * meanG1[i] + n2 * meanG2[i]) / (n1 + n2);
            }
            return -1 * ret;
        }
        // 判別得点
        public double[] score(double[] an, double[][] gx) {
            double[] ret = new double[gx.length];

            for(int i = 0; i < ret.length; i++) {
                ret[i] = an[0] 
                       + an[1] * gx[i][0] 
                       + an[2] * gx[i][1];
            }
            return ret;
        }
    }
}

