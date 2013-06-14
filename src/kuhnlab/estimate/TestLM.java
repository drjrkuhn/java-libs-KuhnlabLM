/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package kuhnlab.estimate;

import com.nr.ch07.GASDEV;
import com.nr.ch15.FGAUSS;
import com.nr.ch15.MRQMIN;
import com.nr.util.COStream;
import com.nr.util.NRMat;
import com.nr.util.NRUtil;
import com.nr.util.NRVec;

/**
 *
 * @author jrkuhn
 */
public class TestLM {

    public static class SumOfGaussians implements EstimateFunction {
        public int nGauss;
        FGAUSS fgauss;
        
        public SumOfGaussians(int numGauss) {
            nGauss = numGauss;
            fgauss = new FGAUSS();
        }
        
        public int getPointDimension() {
            return 1;
        }
        public int getNumCoef() {
            return 3*nGauss;
        }
        public String getCoefName(int index) {
            return null;
        }
        /**
         * y(x; a) is the sum of na/3 Gaussians (15.5.16). The
         * amplitude, center, and width of the Gaussians are stored in
         * consecutive lo cations of a: a[i] = Bk , a[i+1] = Ek ,
         * a[i+2] = Gk , k = 1, ..., na/3. The dimensions of the arrays
         * are a[1..na], dyda[1..na].
         */
        public Estimate getEstimate(double[] point, double[] coef) {
            int numCoef = coef.length;
            double x = point[0];
            double[] ay = {0};
            Estimate est = new Estimate(numCoef);
            fgauss.fgauss(x, coef, ay, est.derivatives);
            est.estimate = ay[0];
            return est;
        }
    }

    //========================================================================
    //========================================================================
    //========================================================================
    
    double a_d[] = {5.0, 2.0, 3.0, 2.0, 5.0, 3.0};
    final double gues_d[] = {4.5, 2.2, 2.8, 2.5, 4.9, 2.8};
    final int NPT = 100,  MA = 6;
    final double SPREAD = 0.001;
    double[] x = NRVec.Double(NPT), y = NRVec.Double(NPT), sig = NRVec.Double(NPT);

    public void initData() {
        GASDEV gasdev = new GASDEV();
        int idum[] = {-911};
        for (int i = 0; i < NPT; i++) {
            x[i] = 0.1 * (i + 1);
            y[i] = 0.0;
            for (int j = 0; j < MA; j += 3) {
                y[i] += a_d[j] * Math.exp(-NRUtil.SQR((x[i] - a_d[j + 1]) / a_d[j + 2]));
            }

            y[i] *= (1.0 + SPREAD * gasdev.gasdev(idum));
            sig[i] = SPREAD * y[i];
        }
    }
    
    public void testLMEstimate() {
        System.out.println("test LMEstimate");

        COStream cout = new COStream(System.out);
        GASDEV gasdev = new GASDEV();

        int idum[] = {-911};
        int i, j;

        SumOfGaussians func = new SumOfGaussians(2);
        LevenbergMarquardtEstimator fit = new LevenbergMarquardtEstimator(func);
        
        // First try a sum of two Gaussians
        for (i = 0; i < NPT; i++) {
            double[] point = new double[] {x[i]};
            fit.addDataPoint(point, y[i], sig[i]);
        }
        fit.setDebugStream(System.out);
        double[] coef = fit.estimate(gues_d);
        double[][] covariance = fit.estimatedCovariance();
        
        cout.print("Results:").endl();
        cout.setw(8).print("a[0]").setw(9).print("a[1]");
        cout.setw(9).print("a[2]").setw(9).print("a[3]");
        cout.setw(9).print("a[4]").setw(9).print("a[5]").endl();
        cout.fixed().setprecision(6);
        for (i = 0; i < 6; i++) {
            cout.setw(9).print(coef[i]);
        }
        cout.endl();
        cout.print("Uncertainties:").endl();
        for (i = 0; i < 6; i++) {
            cout.setw(9).print(Math.sqrt(covariance[i][i]));
        }
        cout.endl().print("Expected results:").endl();
        cout.setw(9).print(5.0).setw(9).print(2.0).setw(9).print(3.0);
        cout.setw(9).print(2.0).setw(9).print(5.0).setw(9).print(3.0).endl();
    }
    
    public void testMRQMin() {
        System.out.println("test MRQMIN");

        COStream cout = new COStream(System.out);
        MRQMIN mrqmin = new MRQMIN();
        FGAUSS fgauss = new FGAUSS();

        int i, j, iter, itst, k, mfit = MA;
        double alamda[] = {0}, chisq[] = {0}, ochisq;
        boolean[] ia = NRVec.Boolean(MA);
        double[] a = NRVec.Double(a_d);
        double[][] covar = NRMat.Double(MA, MA), alpha = NRMat.Double(MA, MA);

        // First try a sum of two Gaussians
        for (i = 0; i < mfit; i++) {
            ia[i] = true;
        }
        for (i = 0; i < MA; i++) {
            a[i] = gues_d[i];
        }
        alamda[0] = -1;
        mrqmin.mrqmin(x, y, sig, a, ia, covar, alpha, chisq, fgauss, alamda);
        k = 1;
        itst = 0;
        for (;;) {
            cout.print("Iteration #").setw(3).print(k);
            cout.print(" test").setw(3).print(itst);
            cout.setw(13).print("chi-squared:").setw(13).print(chisq[0]);
            cout.setw(11).print("alamda[0]:").setw(13).setprecision(10).print(alamda[0]).endl();
            cout.setw(8).print("a[0]").setw(9).print("a[1]");
            cout.setw(9).print("a[2]").setw(9).print("a[3]");
            cout.setw(9).print("a[4]").setw(9).print("a[5]").endl();
            cout.fixed().setprecision(6);
            for (i = 0; i < 6; i++) {
                cout.setw(9).print(a[i]);
            }
            cout . endl();
            k++;

            ochisq = chisq[0];
            mrqmin.mrqmin(x, y, sig, a, ia, covar, alpha, chisq, fgauss, alamda);
            if (Math.abs(ochisq - chisq[0]) < 0.1) {
                itst++;
            } else {
                itst = 0;
            }
            if (itst < 4) {
                continue;
            }
            alamda[0] = 0.0;
            mrqmin.mrqmin(x, y, sig, a, ia, covar, alpha, chisq, fgauss, alamda);
            cout.print("Uncertainties:").endl();
            for (i = 0; i < 6; i++) {
                cout.setw(9).print(Math.sqrt(covar[i][i]));
            }
            cout.endl().print("Expected results:").endl();
            cout.setw(9).print(5.0).setw(9).print(2.0).setw(9).print(3.0);
            cout.setw(9).print(2.0).setw(9).print(5.0).setw(9).print(3.0).endl();
            break;
        }
    }
    
    public static void main(String[] args) {
        TestLM test = new TestLM();
        test.initData();
        test.testMRQMin();
        for (int i=0; i<5; i++) {
            System.out.println("################################################");
        }
        test.testLMEstimate();
    }
    
}
