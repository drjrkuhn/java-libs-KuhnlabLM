/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package kuhnlab.estimate;

import Jama.Matrix;
//import com.nr.util.COStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

/**
 *
 * @author jrkuhn
 */
public class LevenbergMarquardtEstimator {
    protected EstimateFunction function;
    protected ArrayList<Datum> data;
    protected Matrix finalCovariance;

    protected double lambdaStart = 0.001;
    protected double lambdaScalePerStep = 0.1;
    protected double deltaChiSquaredStop = 0.1;
    protected int maxTotalIterations = 1000;
    protected int maxExtraIterations = 3;
    
    protected PrintStream debugStream;

    /**
     * Create a new Estimator to estimate the coefficients of a function.
     * @param function User must provide a class implementing an 
     * ($link EstimateFunction}.
     */
    public LevenbergMarquardtEstimator(EstimateFunction function) {
        data = new ArrayList<Datum>();
        this.function = function;
        debugStream = null;
        finalCovariance = null;
    }

    /**
     * Adds a multidimensional data point to the estimate.
     * 
     * @param point multidimensional point to add. Must have same 
     *              dimensions as {@link EstimateFunction#getPointDimension}
     * @param value value of function at this point
     * @param weight estimated variance of this datapoint. (smaller numbers mean
     *              more weight)
     * @throws java.lang.RuntimeException
     */
    public void addDataPoint(double[] point, double value, double weight) throws RuntimeException {
        if (point.length != function.getPointDimension()) {
            throw new RuntimeException("Data point dimension " + point.length +
                    " does not match function dimension " + function.getPointDimension());
        }
        data.add(new Datum(point, value, weight));
    }
    
    public void setupLambda(double lambdaStart, double lambdaScalePerStep) {
        this.lambdaStart = lambdaStart;
        this.lambdaScalePerStep = lambdaScalePerStep;
    }
    
    public void setupDeltaChiSquaredStopCondition(double deltaChiSquared, int extraIterations) {
        this.deltaChiSquaredStop = deltaChiSquared;
        this.maxExtraIterations = extraIterations;
    }
    
    public void setupMaximumIterations(int maxIterations) {
        this.maxTotalIterations = maxIterations;
    }
    
    public void setDebugStream(PrintStream debug) {
        this.debugStream = debug;
    }

    /**
     * Fit the coefficients to the given datapoints. The user must add 
     * datapoints before calling this function. Use {@link #estimatedCovariance}
     * to get the uncertanties in the final estimate.
     * @param coefGuess Initial guess of the coefficients.
     * @return The final estimated coefficients.
     * @throws java.lang.RuntimeException if the number of coefficients do 
     *      not match, or not enough datapoints have been added, or if the
     *      calculation failes because the matrix of second derivatives
     *      is singular.
     */
    public double[] estimate(double[] coefGuess) throws RuntimeException {
        // Check coefficient dimensions
        int numCoef = coefGuess.length;
        if (numCoef != function.getNumCoef()) {
            throw new RuntimeException("Number of coefficients " + numCoef +
                    " does not match function " + function.getNumCoef());
        }
        // make sure there is enough data to fit
        if (data.size() < numCoef) {
            throw new RuntimeException("There are fewer data points than coefficients");
        }

        // pick a modest value for lambda
        double lambda = lambdaStart;

        Matrix curAlpha = new Matrix(numCoef, numCoef);
        Matrix curBeta = new Matrix(numCoef, 1);
        double[] curCoef = Arrays.copyOf(coefGuess, numCoef);
        
        Matrix trialAlpha = new Matrix(numCoef, numCoef);
        Matrix trialBeta = new Matrix(numCoef, 1);
        double[] trialCoef = new double[numCoef];
        
        Matrix alphaPrime = new Matrix(numCoef, numCoef);
        
        int i, extraIterations = 0, iterations = 1;
        
        // Solve an initial alpha, beta, and chi-sq based on guessed coefficients
        double curChiSq = calcAlphaBetaChi(curCoef, curAlpha, curBeta);
        double lastChiSq = curChiSq;

        do {
            // calculate trial coefficients based on previous alpha and beta
            // and current lambda.

            // (multiply diagonals of alpha by (1+lambda) to form alphaPrime)
            calcAlphaPrime(curAlpha, alphaPrime,lambda);
            
            // (solve the linear equation: AlphaPrime * DeltaCoef = Beta)
            Matrix deltaCoef = alphaPrime.solve(curBeta);
            // apply DeltaCoef to form the coefficients to try.
            for (i = 0; i<numCoef; i++) {
                trialCoef[i] = curCoef[i] + deltaCoef.get(i, 0);
            }
            
            // calculate a new chi-sq (and alpha, beta) based on the trial coefficients
            double trialChiSq = calcAlphaBetaChi(trialCoef, trialAlpha, trialBeta);
            
            if (trialChiSq > curChiSq) {
                if (debugStream != null)
                    debugStream.println("## Bad step. Trying bigger leap.");
                // --Bad step--
                // increase lambda by a factor of 10 and try again with the
                // old (current) alpha and beta
                lambda /= lambdaScalePerStep;
            } else {
                if (debugStream != null)
                    debugStream.println("## On the Golden Path.");
                // --Good step--
                // decrease lambda by a factor of 10 and update coefficients,
                // alpha, and beta to the trial coefficients
                lambda *= lambdaScalePerStep;
                System.arraycopy(trialCoef, 0, curCoef, 0, numCoef);
                Matrices.copyMatrix(trialAlpha, curAlpha);
                Matrices.copyMatrix(trialBeta, curBeta);
            }
            curChiSq = trialChiSq;

            // check to see if the chi-squared value has changed much
            if (Math.abs(curChiSq - lastChiSq) < deltaChiSquaredStop) {
                // we have reached the stop condition, 
                // but allow for a few extra iterations
                extraIterations++;
            } else {
                // we are not at the stop condition
                // clear any extra iterations
                extraIterations = 0;
            }
            
            if (debugStream != null) {
                debugStream.printf("Iteration: %d", iterations);
                if (extraIterations > 0)
                    debugStream.printf(" (%d EXTRA)", extraIterations);
                debugStream.printf(" Chi-Squared: %11.10g%n", curChiSq);
                for (i=0; i<numCoef; i++) {
                    String cname = "coef[" + i + "]";
                    debugStream.printf("%10s ", cname);
                }
                debugStream.println();
                for (i=0; i<numCoef; i++) {
                    debugStream.printf("%10.8g ", curCoef[i]);
                }
                debugStream.println();
            }
            
            lastChiSq = curChiSq;
            iterations++;
        } while (iterations < maxTotalIterations && extraIterations < maxExtraIterations);
        
        if (debugStream != null) {
            debugStream.printf("--- Stopped after %d iterations ---%n", (iterations-1));
        }

        // Compute the covariance matrix as the inverse of the final alpha matrix
        finalCovariance = curAlpha.inverse();
        return curCoef;
    }

    /**
     * Get the covariance of the last estimated coefficients.
     * @return the covariance matrix, or {@code null} if no {@link #estimate}
     *          has not yet been called.
     */
    public double[][] estimatedCovariance() {
        return finalCovariance.getArray();
    }
    
    /**
     * Internal function to calculate alpha, beta, and Chi-Squred
     */
    protected double calcAlphaBetaChi(double[] coef, Matrix alpha, Matrix beta) {
        int row, col, numCoef = coef.length;
        double residual, weightSq;
        double[][] alphaArray = alpha.getArray();
        double[][] betaArray = beta.getArray();
        double chiSq = 0;
        Matrices.clearMatrix(alpha);
        Matrices.clearMatrix(beta);
        // go through each data point and calculate the function and its
        // derivatives with respect to each coefficient at that point.
        for (Datum d : data) {
            EstimateFunction.Estimate est = function.getEstimate(d.point, coef);
            residual = d.value - est.estimate;
            weightSq = d.weight * d.weight;
            chiSq += residual * residual / weightSq;
            // calculate beta (a column vector), the gradient of the function
            // with respect to each coefficient
            for (row = 0; row < numCoef; row++) {
                betaArray[row][0] += residual * est.derivatives[row] / weightSq;
            }
            // calculate alpha, the [approximate] Hessian matrix
            for (row = 0; row < numCoef; row++) {
                for (col = 0; col < numCoef; col++) {
                    alphaArray[row][col] += est.derivatives[row] 
                            * est.derivatives[col] / weightSq;
                }
            }
        }
        return chiSq;
    }
    
    /**
     * Internal function to calculate alphaPrime from alpha and lambda
     */
    protected void calcAlphaPrime(Matrix alpha, Matrix alphaPrime, double lambda) {
        Matrices.copyMatrix(alpha, alphaPrime);
        double[][] alphaPrimeArray = alphaPrime.getArray();
        int numCoef = alpha.getRowDimension();
        double onePlusLambda = 1.0 + lambda;
        // Adjust the diagonals of a by (1 + lambda)
        for (int diag=0; diag<numCoef; diag++) {
            alphaPrimeArray[diag][diag] *= onePlusLambda;
        }
    }

    /**
     * Internal class used to store data.
     */
    protected static class Datum {
        protected double[] point;
        protected double value;
        protected double weight;

        protected Datum(double[] point, double value, double weight) {
            this.point = point;
            this.value = value;
            this.weight = weight;
        }
    }
}
