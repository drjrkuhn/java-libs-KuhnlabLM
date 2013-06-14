/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package kuhnlab.estimate;

/** 
 * Used to implement a function, f(<b>X</b>;<b>A</b>). Where the variables <B>X</b> and
 * function coefficients <b>A</b> are both vectors.
 * <p>
 * See {@link LevenbergMarquardtEstimator}.
 *
 * @author jrkuhn
 */
public interface EstimateFunction {
    /**
     * Return value for getEstimate.
     */
    public static class Estimate {
        /** Estimated value of {@code f} at the requested point <b>X</b>. */
        public double estimate;
        /** Derivatives of f with respect to each parameter, <b>A</b>. */
        public double[] derivatives;
        
        public Estimate(int numCoeff) {
            derivatives = new double[numCoeff];
        }
        
        public String toString() {
            String str = "estimate = " + estimate + "\n";
            for (int i=0; i<derivatives.length; i++) {
                str += "derivative["+i+"] = " + derivatives[i] + "\n";
            }
            return str;
        }
    }
    /**
     * Get number of "x" values. For example: 1 for f(x;A), 2 for f(x,y;A), etc.
     * @return number of dimensions of a point
     */
    public int getPointDimension();

    /**
     * Get number of coefficients a that determine this function f(x;a).
     * @return number of coefficients in in f
     */
    public int getNumCoef();
    
    /**
     * Get the name of a coefficent.
     * @param index index into array of coefficents
     * @return the name of the coefficient.
     */
    public String getCoefName(int index);
    
    /**
     * Get the estimated function value and coefficients at a point.
     * @param point     point to evalute function at
     * @param coef      coefficients to use to evaluate function
     * @return the estimated value and derivatives at the data point
     */
    public Estimate getEstimate(double[] point, double coef[]);
}
