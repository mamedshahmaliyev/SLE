package sle_test;

import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.concurrent.TimeUnit;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class SLE_TEST {
    private double lyambda = 30;
    private double mu1 = 15;
    private double mu2 = 3;
    private double sigma1 = 0.3;
    private double sigma2 = 0.7;
    private double fi1 = 0.6;
    private double fi2 = 0.4;
    private double gamma = 2;
    private double nu = 1;
    private double tao = 0.5;
    private int s = 10;
    private int S = 50;
    private int N = 80;
    
    private long getDateDiff(Date startDate, Date endDate, TimeUnit timeUnit) {
        long diffInMillies = endDate.getTime() - startDate.getTime();
        return timeUnit.convert(diffInMillies,TimeUnit.MILLISECONDS);
    }
    private String getDateDiffStr(Date startDate, Date endDate) {
        return getDateDiff(startDate, endDate, TimeUnit.HOURS) +" hours or "+getDateDiff(startDate, endDate, TimeUnit.MINUTES) +" minutes or "+getDateDiff(startDate, endDate, TimeUnit.SECONDS) +" seconds or "+getDateDiff(startDate, endDate, TimeUnit.MILLISECONDS) +" milliseconds ";
    }
    
    private void solve(double[][] A, double[] B){
        RealMatrix n = new Array2DRowRealMatrix(A);
        RealVector constants = new ArrayRealVector(B, false);
        Date startDate = new Date();
        String format = "yyyy.MM.dd hh:mm:sss.S";
        System.out.println((new SimpleDateFormat (format)).format(startDate));
        DecompositionSolver solver = new LUDecomposition(n).getSolver();
        solver.solve(constants);
        Date endDate = new Date();
        System.out.println((new SimpleDateFormat (format)).format(endDate));
        System.out.println("S * N = " +S + " * " + N +" = " + (S * N)+": "+getDateDiffStr(startDate, endDate));
        //return solver.solve(constants).toArray();
    }
    
    private double getTransitionMatrixElement(int m1, int n1, int m2, int n2){
        //VSO policy
        if (m1 > s){
            if (m2 == m1 && n2 == n1 + 1) return lyambda;
            else if (m2 == m1 && n2 == n1 - 1) return mu1 * sigma1;
            else if (m2 == m1 - 1 && n2 == n1 - 1) return mu2 * sigma2;
            else if (m2 == m1 - 1 && n2 == n1 && n2 == 0) return m1 * gamma;
            else if (m2 == m1 - 1 && n2 == n1 && n2 > 0) return (m1 - 1) * gamma;
        } else
        if (m1 <= s && m1 > 0){
            if (m2 == m1 && n2 == n1 + 1) return lyambda;
            else if (m2 == m1 && n2 == n1 - 1) return mu1 * sigma1;
            else if (m2 == m1 - 1 && n2 == n1 - 1) return mu2 * sigma2;
            else if (m2 == m1 - 1 && n2 == n1 && n2 == 0) return m1 * gamma;
            else if (m2 == m1 - 1 && n2 == n1 && n2 > 0) return (m1 - 1) * gamma;
            else if (m2 == S && n2 == n1) return nu;
        } else
        if (m1 == 0){
            if (m2 == 0 && n2 == n1 + 1) return lyambda * fi1;
            else if (m2 == 0 && n2 == n1 - 1) return n1 * tao;
            else if (m2 == S && n2 == n1) return nu;
        }
        return 0;
    }
    public static void main(String[] args) {
        SLE_TEST sle = new SLE_TEST();
        
        int S_MIN = 50;
        int S_MAX = 100;
        int N_MIN = 50;
        int N_MAX = 100;
        
        for (int m = S_MIN; m <= S_MAX; m = m + 20) {
            for (int n = N_MIN; n <= N_MAX; n = n + 20) {
                //set N, S and s
                sle.S = m;
                sle.s = sle.S / 2 - 1;
                sle.N = n;
                //generate transition matrix (transposed and last row replaced with 1) 
                int last_row =  sle.S*(sle.N+1)+sle.N;
                double[][] TM = new double[(sle.S+1) * (sle.N+1)][(sle.S+1) * (sle.N+1)];
                for (int m1 = 0; m1 <= sle.S; m1++) {
                    for (int n1 = 0; n1 <= sle.N; n1++) {
                        int i = m1*(sle.N+1)+n1;
                        double row_sum = 0;
                        for (int m2 = 0; m2 <= sle.S; m2++) {
                            for (int n2 = 0; n2 <= sle.N; n2++) {
                                int j = m2*(sle.N+1)+n2;
                                TM[j][i] = sle.getTransitionMatrixElement(m1, n1, m2, n2);
                                row_sum += TM[j][i]; 
                                if (j == last_row) TM[j][i] = 1;
                            }
                        }
                        //set diagonal elements as negative row sum
                        TM[i][i] = -1 * row_sum;
                    }
                }
                
                double[] b = new double[TM.length];
                b[b.length-1] = 1;
                
                //find steady state probabilities as solution of corresponding SLE
                sle.solve(TM, b);
            }
        }
    }
    
}
