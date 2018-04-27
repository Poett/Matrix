import java.util.Random;

public class Matrix {
	
	public static int operationCount;
	private double[][] L;
	private double[][] U;
	
	 // return a random m-by-n matrix with values between 0 and 1
    public static double[][] random(int m, int n, int bound) {
        double[][] a = new double[m][n];
        Random rand = new Random();
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                a[i][j] = rand.nextInt(bound);
        return a;
    }
    
    public static double[] random(int m) {
        double[] a = new double[m];
        Random rand = new Random();
        for (int i = 0; i < m; i++)
        	a[i] = rand.nextDouble();
        return a;
    }

    // return n-by-n identity matrix I
    public static double[][] identity(int n) {
        double[][] a = new double[n][n];
        for (int i = 0; i < n; i++)
            a[i][i] = 1;
        return a;
    }

    // return x^T y
    public static double dot(double[] x, double[] y) {
        if (x.length != y.length) throw new RuntimeException("Illegal vector dimensions.");
        double sum = 0.0;
        for (int i = 0; i < x.length; i++)
            sum += x[i] * y[i];
        return sum;
    }

    // return B = A^T
    public static double[][] transpose(double[][] a) {
        int m = a.length;
        int n = a[0].length;
        double[][] b = new double[n][m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                b[j][i] = a[i][j];
        return b;
    }

    // return c = a + b
    public static double[][] add(double[][] a, double[][] b) {
        int m = a.length;
        int n = a[0].length;
        double[][] c = new double[m][n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                c[i][j] = a[i][j] + b[i][j];
        return c;
    }

    // return c = a - b
    public static double[][] subtract(double[][] a, double[][] b) {
        int m = a.length;
        int n = a[0].length;
        double[][] c = new double[m][n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                c[i][j] = a[i][j] - b[i][j];
        return c;
    }

    // return c = a * b
    public static double[][] multiply(double[][] a, double[][] b) {
        int m1 = a.length;
        int n1 = a[0].length;
        int m2 = b.length;
        int n2 = b[0].length;
        if (n1 != m2) throw new RuntimeException("Illegal matrix dimensions.");
        double[][] c = new double[m1][n2];
        for (int i = 0; i < m1; i++)
            for (int j = 0; j < n2; j++)
                for (int k = 0; k < n1; k++)
                    c[i][j] += a[i][k] * b[k][j];
        return c;
    }

    // matrix-vector multiplication (y = A * x)
    public static double[] multiply(double[][] a, double[] x) {
        int m = a.length;
        int n = a[0].length;
        if (x.length != n) throw new RuntimeException("Illegal matrix dimensions.");
        double[] y = new double[m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                y[i] += a[i][j] * x[j];
        return y;
    }


    // vector-matrix multiplication (y = x^T A)
    public static double[] multiply(double[] x, double[][] a) {
        int m = a.length;
        int n = a[0].length;
        if (x.length != m) throw new RuntimeException("Illegal matrix dimensions.");
        double[] y = new double[n];
        for (int j = 0; j < n; j++)
            for (int i = 0; i < m; i++)
                y[j] += a[i][j] * x[i];
        return y;
    }

    
    public static String toString(double[][] a) 
    {
    	String toReturn = "";
    	
    	for(int i = 0; i < a[0].length; i++) 
    	{
    		for(int j = 0; j < a.length; j++) 
    		{
    			toReturn += a[i][j] + "\t";
    		}
    		toReturn += "\n";
    	}
    	
    	return toReturn;
    }
    
    public static String toString(double[] a) 
    {
    	String toReturn = "";
    	
    	for(int i = 0; i < a.length; i++) 
    	{
    		toReturn += a[i] + "\n";
    	}
    	
    	return toReturn;
    }
    
    public static double[][] decomposeLU(double[][] A)
    {
    	

    	int n = A.length;
    	double[][] LU = new double[n][n];

    	
        for (int p = 0; p < n; p++) {
            //Gauss elimination to echelon
            for (int i = p + 1; i < n; i++) {
                double alpha = A[i][p] / A[p][p];
                if(alpha != 0) LU[i][p] = alpha;
                for (int j = p; j < n; j++) {
                    A[i][j] -= alpha * A[p][j];
                    operationCount++;
                }
            }
        }
        //copy upper triangular to LU except 1's
        int j = 0;
        int k;
        while(j < n) 
        {
        	k = j;
        	while(k < n) 
        	{
        		LU[j][k] = A[j][k];
        		k++;
        	}
        	j++;
        }
    	return LU;
    }
    
    public static double[] solveLU(double[][] A, double[] b) 
    {
    	int n = A.length;
    	
    	double[][] U = new double[n][n];
    	double[][] L = Matrix.identity(n);
    	
    	//split the stored LU comp
    	for(int i = 0; i < n; i++) 
    	{
    		for(int j = 0; j < n; j++) 
    		{
    			if(i <= j) {U[i][j] = A[i][j];} //if we are in diagonal or above triangle
    			else {
    				L[i][j] = A[i][j];
    				U[i][j] = 0;
    				}
    		}
    	}
    	
    	//Forward substitution on L
    	//Backward Substitution on U
    	return Matrix.backwardSub(U, Matrix.forwardSub(L, b));
    }
    
    public static double[] solveGauss(double[][] A, double[] b) 
    {
    	int n = b.length;

        for (int p = 0; p < n; p++) {
            // pivot within A and b
            for (int i = p + 1; i < n; i++) {
                double alpha = A[i][p] / A[p][p];
                b[i] -= alpha * b[p];
                operationCount++;
                for (int j = p; j < n; j++) {
                    A[i][j] -= alpha * A[p][j];
                    operationCount++;
                }
            }
        }
        // back substitution
        return Matrix.backwardSub(A, b);

    }
    
    public static double[] backwardSub(double[][] A, double[] b) 
    {
    	int n = A.length;
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < n; j++) {
                sum += A[i][j] * x[j];
                operationCount++;
            }
            x[i] = (b[i] - sum) / A[i][i];
            operationCount++;
        }
        return x;
    }
    
    public static double[] forwardSub(double[][] A, double[] b) 
    {
    	int n = A.length;
        double[] x = new double[n];
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                sum += A[i][j] * x[j];
                operationCount++;
            }
            x[i] = (b[i] - sum) / A[i][i];
            operationCount++;
        }
        return x;
    }
    
    public static void print(double[][] a) 
    {
    	System.out.println(toString(a));
    }
    
    public static void print(double[] a) 
    {
    	System.out.println(toString(a));
    }
}

