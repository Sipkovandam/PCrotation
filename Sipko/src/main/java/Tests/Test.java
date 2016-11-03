package Tests;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.stream.IntStream;

import org.junit.rules.TemporaryFolder;

import PCA.Matrix;
import PCA.MatrixStruct;
import PCA.RLog;
import PCA_testcases.CompareFiles;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

public class Test {
	
	public static void main(String[] args) throws IOException
	{	
		double testval= 1.7;
		
		double[] testVals = new double[]{1,2,3,4,5};
		IntStream.range(0, testVals.length).forEach(x -> System.out.println(testVals[x]));


//		String geoFN= "E:/Groningen/Scripts/Tests/Rlog.java/DESeqNorm/Samples.DESeqNorm.txt.gz";
//		CompareFiles.compare(geoFN,"TestCaseFiles/Result_Samples.txt",true);
//		

	}
	public void print()
	{
		
		
	}
}
