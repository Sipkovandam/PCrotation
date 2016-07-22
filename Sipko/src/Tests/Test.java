package Tests;

import PCA.Matrix;
import pca.MatrixStruct;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

public class Test {
	
	public static void main(String[] args)
	{
		Matrix test = new Matrix("E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/TESTexpression_correlation.txt");
		MatrixStruct test2 = new MatrixStruct("E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/TESTexpression_correlation.txt");
		
		
		MatrixStruct test4 = new MatrixStruct("E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/TESTexpression_correlation.txt.gz");
		Matrix test3 = new Matrix("E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/TESTexpression_correlation.txt.gz");
		System.out.println(test3.values[2][3]);
	}

}
