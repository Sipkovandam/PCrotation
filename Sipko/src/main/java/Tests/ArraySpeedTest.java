package Tests;

import org.ujmp.core.DenseMatrix;

//import no.uib.cipr.matrix.DenseMatrix;

import org.ujmp.core.Matrix;
import org.ujmp.core.SparseMatrix;

//import org.ujmp.core.doublematrix.SparseDoubleMatrix2D;
//import org.ujmp.mtj.MTJDenseDoubleMatrix2D;
//
import no.uib.cipr.matrix.EVD;
//
//import org.ujmp.core.DenseMatrix;
//
//import org.ujmp.core.Matrix;
//import org.ujmp.core.bigdecimalmatrix.BigDecimalMatrix;
//import org.ujmp.core.doublematrix.DenseDoubleMatrix2D;
import org.ujmp.core.doublematrix.calculation.general.decomposition.SVD;
//import org.ujmp.core.util.UJMPSettings;
//import org.ujmp.core.util.matrices.MatrixLibraries;
//import org.ujmp.mtj.MTJDenseDoubleMatrix2D;


//import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.EVD;
//import no.uib.cipr.matrix.NotConvergedException;


public class ArraySpeedTest
{
	//A test to compare speeds. Only works with squared array as I was to lazy to figure out which dimensions each matrix should have
	private class MatrixStandard
	{
		double[][] values = null;
		MatrixStandard(int rows, int cols)
		{
			values = new double[rows][cols];
		}
		public int numRows()
		{
			return this.values.length;
		}
		public int numColumns()
		{
			return this.values[0].length;
		}
	}
	private class MatrixStruct
	{
		double[] values = null;
		int numRows = 0;
		int numColumns = 0;
		
		MatrixStruct(int rows, int cols)
		{
			values = new double[rows*cols];
			numRows=rows;
			numColumns=cols;
		}
		public double getVal(int row, int col)
		{
			return this.values[numColumns*row+col];

		}
		public void setVal(int row, int col, double value)
		{
			this.values[numColumns*row+col]=value;
		}
		
		public int numRows()
		{
			return this.numRows;
		}
		public int numColumns()
		{
			return this.numColumns;
		}
	}
	
	public static void main(String[] args) throws Exception
	{
		ArraySpeedTest arraySpeedTest = new ArraySpeedTest();
		int rows=32222;
		int cols=52222;
		System.out.println("rows = " + rows + " cols=" + cols);
//		MatrixStruct matrixStruct = arraySpeedTest.new MatrixStruct(rows, cols);
//		MatrixStandard matrix = arraySpeedTest.new MatrixStandard(rows, cols);
		System.out.println("Bla1234");
		//no.uib.cipr.matrix.DenseMatrix denseMatrix = new no.uib.cipr.matrix.DenseMatrix(rows, cols);
		System.out.println("Blabla");
		
//		SparseMatrix m1 = SparseMatrix.Factory.zeros(rows,cols);
//		
//		for(int r= 0; r < rows; r++)
//		{
//			if(r%100==0)
//				System.out.println("row =" + r + "/"+ rows);
//			for(int c = 0; c < cols; c++)
//			{
//				double val = Math.random();
//				m1.setAsDouble(val, r,c);
//			}
//		}
		System.out.println("Done!");
		Matrix.Factory.availableProcessors();

		Matrix rand = Matrix.Factory.zeros(rows,cols);
		System.out.println("Blablabla");
		
		for(int r= 0; r < rows; r++)
			for(int c = 0; c < cols; c++)
			{
				double val = Math.random();
//				matrix.values[r][c]=val;
//				matrixStruct.setVal(r,c,val);
				rand.setAsDouble(val, r,c);
//				denseMatrix.set(r, c, val);
			}
//		System.out.println("1()=" + UJMPSettings.getInstance().isUseEJML());
//		no.uib.cipr.matrix.DenseMatrix m = null;
//		System.out.println("2()=" + UJMPSettings.getInstance().isUseEJML());
//		if (rand instanceof MTJDenseDoubleMatrix2D) {
//			m = ((MTJDenseDoubleMatrix2D) rand).getWrappedObject();
//		} else {
//			m = new MTJDenseDoubleMatrix2D(rand).getWrappedObject();
//		}
//		System.out.println("3()=" + UJMPSettings.getInstance().isUseEJML());
////		no.uib.cipr.matrix.SVD svd = no.uib.cipr.matrix.SVD.factorize(m);
//		System.out.println("EJML()=" + UJMPSettings.getInstance().isUseEJML());
//		System.out.println("MTJ()=" + UJMPSettings.getInstance().isUseMTJ());

		
//		MatrixStandard Z = arraySpeedTest.new MatrixStandard(matrix.values.length, matrix.values.length);
		long start0 = System.nanoTime();
		//Matrix res2 = rand.mtimes(rand);
		Matrix res2 = SVD.INSTANCE.calc(rand)[0];
//		Matrix res2 = rand.eig()[0];
		long end0 = System.nanoTime();
		System.out.println("runtime=" + (end0-start0));
		System.out.println("runtime=" + (end0-start0)/1000+ "nanosec");
		System.out.println("runtime=" + (end0-start0)/1000000 + "ms");
		System.out.println("runtime=" + (end0-start0)/1000000000 + "sec");
		
//		System.out.println("\nstd");
//		MatrixStandard C = arraySpeedTest.new MatrixStandard(matrix.values.length, matrix.values.length);
//		long start = System.nanoTime();
//		multAdd(matrix,matrix,C);
//		long end = System.nanoTime();
//		System.out.println("runtime=" + (end-start));
//		System.out.println("runtime=" + (end-start)/1000+ "nanosec");
//		System.out.println("runtime=" + (end-start)/1000000 + "ms");
//		System.out.println("runtime=" + (end-start)/1000000000 + "sec");
		
//		System.out.println("\n1 array");
//		MatrixStruct D = arraySpeedTest.new MatrixStruct(matrix.values.length, matrix.values.length);
//		long start2 = System.nanoTime();
//		multAdd2(matrixStruct,matrixStruct,D);
//		long end2 = System.nanoTime();
//		System.out.println("runtime=" + (end2-start2));
//		System.out.println("runtime=" + (end2-start2)/1000+ "nanosec");
//		System.out.println("runtime=" + (end2-start2)/1000000 + "ms");
//		System.out.println("runtime=" + (end2-start2)/1000000000 + "sec");
	
		
		
		
		
		
		
		
		EVD evd = new EVD(rows, false, true);
        
		System.out.println("\n DenseMatrix");
		//DenseMatrix E = new DenseMatrix(matrix.values.length, matrix.values.length);
		long start3 = System.nanoTime();
		//denseMatrix.mult(denseMatrix, E);
	//	no.uib.cipr.matrix.DenseMatrix res = evd.factor(denseMatrix).getRightEigenvectors();
		long end3 = System.nanoTime();
		System.out.println("runtime=" + (end3-start3));
		System.out.println("runtime=" + (end3-start3)/1000+ "nanosec");
		System.out.println("runtime=" + (end3-start3)/1000000 + "ms");
		System.out.println("runtime=" + (end3-start3)/1000000000 + "sec");

		
		
		
		
		
		
		
		
		
		
		
		
		
		//		
//		System.out.println("\n DenseMatrix2");
//		DenseMatrix F = new DenseMatrix(matrix.values.length, matrix.values.length);
//		long start4 = System.nanoTime();
//		multAdd3(denseMatrix,denseMatrix,F);
//		long end4 = System.nanoTime();
//		System.out.println("runtime=" + (end4-start4));
//		System.out.println("runtime=" + (end4-start4)/1000+ "nanosec");
//		System.out.println("runtime=" + (end4-start4)/1000000 + "ms");
//		System.out.println("runtime=" + (end4-start4)/1000000000 + "sec");
		
//		String output = "";
//		String input = "";
//		for(int r = 0; r< 5 ; r++)//res.numRows(); r++)
//		{
//			String line = Double.toString(res.get(r,0));
//			String inLine = Double.toString(res2.getAsDouble(r,0));
//			for(int c = 0; c < 5;c++)//res.numColumns(); c++)
//			{
//				line += "\t" + res.get(r,c);
//				inLine += "\t" + res2.getAsDouble(r,c);
//				//System.out.println("newVS old=\t" + res.get(r,c) + "\t" + res3[r]);
//				//System.out.println("newVS old=\t" + res2.getAsDouble(r,c) + "\t" + E.get(r,c));
//			}
//			output+=line+"\n";
//			input+=inLine+"\n";
//		
//		}
//		System.out.println("input\n"+input);
//		System.out.println("output\n"+output);
	}
	public static MatrixStandard multAdd(MatrixStandard A, MatrixStandard B, MatrixStandard C) {
		
		int numRows = A.numRows();
		int numColumns = A.numColumns();
        for (int i = 0; i < numRows; ++i)
            for (int j = 0; j < C.numColumns(); ++j) {
                double dot = 0;
                for (int k = 0; k < numColumns; ++k)
                    dot += A.values[i][k] * B.values[k][j];
                C.values[i][j]+=dot;
            }

        return C;
    }
	
	public static MatrixStruct multAdd2(MatrixStruct A, MatrixStruct B, MatrixStruct C) {
		
		int numRows = A.numRows();
		int numColumns = A.numColumns();
        for (int i = 0; i < numRows; ++i)
            for (int j = 0; j < C.numColumns(); ++j) {
                double dot = 0;
                for (int k = 0; k < numColumns; ++k)
                    dot += A.getVal(i,k) * B.getVal(k,j);
                C.setVal(i,j, dot);
            }
        return C;
    }
//	public static DenseMatrix multAdd3(DenseMatrix A, DenseMatrix B, DenseMatrix C) {
//		
//		int numRows = A.numRows();
//		int numColumns = A.numColumns();
//        for (int i = 0; i < numRows; ++i)
//            for (int j = 0; j < C.numColumns(); ++j) {
//                double dot = 0;
//                for (int k = 0; k < numColumns; ++k)
//                    dot += A.get(i,k) * B.get(k,j);
//                C.set(i,j, dot);
//            }
//        return C;
//    }
}
