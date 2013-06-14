/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package kuhnlab.estimate;

import Jama.Matrix;
import java.util.Arrays;

/**
 * Utility class for JAMA Matrix library. This class adds some extra functionality
 * missing from the JAMA Matrix library.
 * 
 * @author jrkuhn
 */
public class Matrices {
    
    /**
     * Clear (set to zero) the contents of a JAMA Matrix.
     * @param A the Matrix to clear
     */
    public static void clearMatrix(Matrix A) {
        double[][] array = A.getArray();
        int m = A.getRowDimension();
        for (int row = 0; row < m; row++) {
            Arrays.fill(array[row], 0.0);
        }
    }
    
    /**
     * Copy data from a source matrix to a destination Matrix.
     * @param Src source Matrix to copy.
     * @param Dest destination Matrix to fill.
     */
    public static void copyMatrix(Matrix Src, Matrix Dest) {
        double[][] sarray = Src.getArray();
        double[][] darray = Dest.getArray();
        int m = Src.getRowDimension();
        int n = Src.getColumnDimension();
        for (int row=0; row<m; row++) {
            System.arraycopy(sarray[row], 0, darray[row], 0, n);
        }
    }


}
