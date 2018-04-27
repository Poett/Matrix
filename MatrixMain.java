import java.util.ArrayList;

public class MatrixMain {

    public static void main(String[] args) {
    	
    	double[][] g = {
    			{1, 3, 4},
    			{2, 1, -1},
    			{4, 2, -1}
    	};
    	
    	Matrix.print(Matrix.decomposeLU(g));
    	
    	
    	double[][] R = Matrix.random(15, 15, 10);
    	Matrix.print(R);
    	
    	
    	double[][] A = {
    			{1,5,2,-3,-4},
    			{2,1,-2,-1,0},
    			{0,2,5,3,-2,},
    			{0,8,-3,2,-1},
    			{1,0,3,-2,2}
    	};
        double[][] A2 = {
    			{1,5,2,-3,-4},
    			{2,1,-2,-1,0},
    			{0,2,5,3,-2,},
    			{0,8,-3,2,-1},
    			{1,0,3,-2,2}
    	};
        
        A = R.clone();
        A2 = R.clone();
        

        
        ArrayList<double[]> b = new ArrayList<>();
        ArrayList<double[]> c = new ArrayList<>();
        for(int i = 0; i < 1024; i++) 
        {
        	double[] random = Matrix.random(15);
        	b.add(random);
        	c.add(random);
        }
        
        
        
        //Gaussian Testing
        Matrix.operationCount = 0;
        long start = System.currentTimeMillis();
        for(double[] solution : b) 
        {
        	Matrix.solveGauss(A, solution);
        }
        long end = System.currentTimeMillis();
        System.out.println("--------Gaussian---------");
        System.out.println("Operations for 1024 linear system solutions of a 5x5 matrix: " + Matrix.operationCount + " in " + (end - start) + "miliseconds");
        
   
        //LU Decomposition

        Matrix.operationCount = 0;
        double[][] LU = Matrix.decomposeLU(A2);
        
        start = System.currentTimeMillis();
        for(double[] solution : c) 
        {
        	Matrix.solveLU(LU, solution);
        }
        end = System.currentTimeMillis();
        
        System.out.println("--------LU Factorization---------");
        System.out.println("Operations for 1024 linear system solutions of a 5x5 matrix: " + Matrix.operationCount + " in " + (end - start) + "miliseconds");
        

        
        
        
    }

}
