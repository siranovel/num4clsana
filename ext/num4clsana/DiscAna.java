import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.stat.correlation.Covariance;
import org.apache.commons.math3.stat.StatUtils;

import java.util.Map;
import java.util.HashMap;
public class DiscAna {
    private static DiscAna da = new DiscAna();
    public static DiscAna getInstance() { return da;}
    public Map<String, double[]> score(double[][] xa, double[][] xb){
        LineDisc da0 = new LineDisc();

        return da0.score(xa, xb);
    }
    public Map<String, double[]> score2(double[][] xa, double[][] xb) {
        Mahalanobis da1 = new Mahalanobis();

        return da1.score(xa, xb);
    }
    public Map<String, Double> validity(Map<String, double[]> sc, double[][] xa, double[][] xb) {
        Validity va = new Validity();

        return va.validity(sc, xa, xb);
    }
    /*********************************/
    /* interface define              */
    /*********************************/
    /*********************************/
    /* class define                  */
    /*********************************/
    // 線形型判別分析
    private class LineDisc {
        private int n1 = 0;
        private int n2 = 0;
        // 判別得点
        public Map<String, double[]> score(double[][] xa, double[][] xb) {
            Map<String, double[]> retMap = new HashMap<String, double[]>();
            double[] an = line_disc(xa, xb);  

            retMap.put("G1", calcScore(an, xa));
            retMap.put("G2", calcScore(an, xb));
            return retMap;
        }
        // 線形型判別関数の係数
        private double[] line_disc(double[][] xa, double[][] xb){
            // 各グループの個数
            RealMatrix matrixXa = MatrixUtils.createRealMatrix(xa);
            n1 = xa.length;
            RealMatrix matrixXb = MatrixUtils.createRealMatrix(xb);
            n2 = xb.length;

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
        // 平均を計算
        private double[] calcMean(RealMatrix matrixX) {
            double[] m = new double[2];

            for(int i = 0; i < 2; i++) {
                double[] dt = matrixX.getColumn(i);

                m[i] = StatUtils.mean(dt);
            }
            return m;
        }
        // 分散共分散行列を計算
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
        private double[] calcScore(double[] an, double[][] gx) {
            double[] ret = new double[gx.length];

            for(int i = 0; i < ret.length; i++) {
                ret[i] = an[0] 
                       + an[1] * gx[i][0] 
                       + an[2] * gx[i][1];
            }
            return ret;
        }
    }
    // マハラノビスの距離
    private class Mahalanobis {
        // 判別得点
        public Map<String, double[]> score(double[][] xa, double[][] xb) {
            Map<String, double[]> retMap = new HashMap<String, double[]>();
            RealMatrix matrixXa = MatrixUtils.createRealMatrix(xa);
            RealMatrix matrixXb = MatrixUtils.createRealMatrix(xb);
            // 各グループの平均を計算
            double[] meanG1 = calcMean(matrixXa);
            double[] meanG2 = calcMean(matrixXb);
            // 各グループの分散共分散行列を求め、逆行列を計算
            RealMatrix invV1 = MatrixUtils.inverse(calcCorrMatrix(matrixXa));
            RealMatrix invV2 = MatrixUtils.inverse(calcCorrMatrix(matrixXb));

            // 判別得点を計算
            retMap.put("G1", 
                calcScore(
                    xa,
                    meanG1, meanG2,
                    invV1, invV2)
            );
            retMap.put("G2", 
                calcScore(
                    xb,
                    meanG1, meanG2,
                    invV1, invV2)
            );
            return retMap;
        }
        // 平均を計算
        private double[] calcMean(RealMatrix matrixX) {
            double[] m = new double[2];

            for(int i = 0; i < 2; i++) {
                double[] dt = matrixX.getColumn(i);

                m[i] = StatUtils.mean(dt);
            }
            return m;
        }
        // 分散共分散行列を計算
        private RealMatrix calcCorrMatrix(RealMatrix matrixX) {
            Covariance corel = new Covariance();
            double[][] corr = new double[2][2];

            for(int i =0; i < 2; i++) {
                for(int j = 0; j < 2; j++) {
                    double[] xArray = matrixX.getColumn(i);
                    double[] yArray = matrixX.getColumn(j);
                    
                    corr[i][j] = corel.covariance(xArray, yArray);
                }
            }
            return MatrixUtils.createRealMatrix(corr);
        }
        private double[] calcScore(
            double[][] xi,
            double[] meanG1, double[] meanG2,
            RealMatrix invV1, RealMatrix invV2) {
            double[] ret = new double[xi.length];

            for(int i = 0; i < xi.length; i++) {
                double d1 = calcD(xi[i], meanG1,  invV1);              
                double d2 = calcD(xi[i], meanG2,  invV2); 

                ret[i] = Math.sqrt(d2) - Math.sqrt(d1);             
            }
            return ret;
        }
        private double calcD(double[] x, double[] mean, RealMatrix InvV) {
            RealMatrix xVec = MatrixUtils.createRealMatrix(
                new double[][] {
                    {x[0] - mean[0]},
                    {x[1] - mean[1]}
               }
            );
            RealMatrix xt = xVec.transpose();
            // tx*(v-1)*x
            RealMatrix dMatrix = xt.multiply(InvV.multiply(xVec));

            return new LUDecomposition(dMatrix).getDeterminant();
        }
    }
    // 
    private class Validity {
        private int act_0 = 0;         // 実測値 0
        private int act_1 = 0;         // 実測値 1
        private int predict_0 = 0;
        private int predict_1 = 0;
        private double tp = 0.0;       // 実測値=0 予測値 = 1
        private double fn = 0.0;       // 実測値=0 予測値 = 0
        private double fp = 0.0;       // 実測値=1 予測値 = 1
        private double tn = 0.0;       // 実測値=1 予測値 = 0
        public Map<String, Double> validity(Map<String, double[]> sc, double[][] xa, double[][] xb) {
            resClass(sc, xa, xb);      // 分類結果
            return preciIndex();       // 分類問題の精度指標
        }
        // 分類結果
        private void resClass(Map<String, double[]> sc, double[][] xa, double[][] xb) {
            act_0 = xa.length;
            act_1 = xb.length;
            // 実測値=0
            calcPredict(sc.get("G1"));
            tp = (double)predict_0 / (double)act_0;
            fn = (double)predict_1 / (double)act_0;
            // 実測値=1
            calcPredict(sc.get("G2"));
            fp = (double)predict_0 / (double)act_1;
            tn = (double)predict_1 / (double)act_1;
        }
        // 予測値
        private void calcPredict(double[] gn) {
            predict_0 = 0;
            predict_1 = 0;

            for(double v : gn) {
                if (v > 0)  { predict_0++;}
                if (v < 0)  { predict_1++;}
            }
        }
        // 分類問題の精度指標
        private Map<String, Double> preciIndex() {
            Map<String, Double> retMap = new HashMap<String, Double>();

            // 精度・正確度
            retMap.put("accuracy", (tp+tn)/(tp+fp+tn+fn));
            // 適合度
            retMap.put("precision", tp/(tp+fp));
            // 再現度
            retMap.put("recall", tp/(tp+fn));
            // 感度=再現率
            retMap.put("sensitivity", tp/(tp+fn));
            // 特異度
            retMap.put("specificity", tn/(fp+tn));
            // 陽性敵中率=適合度
            retMap.put("ppv", tp/(tp+fp));
            // 陰性的中率
            retMap.put("npv", tn/(tn+fn));
            // 真陽性率=再現率・感度
            retMap.put("tpr", tp/(tp+fn));
            // 偽陽性率
            retMap.put("fpr", fp/(fp+tn));
            // 真陰性率
            retMap.put("tnr", tn/(fp+tn));
            // 偽陰性率
            retMap.put("fnr", fn/(tp+fn));
            return retMap;
        }
    }

}

