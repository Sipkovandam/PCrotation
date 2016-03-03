package PCA;

import java.io.IOException;
import java.util.Arrays;
import java.util.stream.Collectors;

import pca.MatrixStruct;

public class test 
{
	
	public static void main(String[] args) throws IOException
	{
		String testFN = "E:/Groningen/Data/LifeLines/Phenotypes/Kallisto mapping counts/AC1C40ACXX-1-18_220/test.txt";
		MatrixStruct test = new MatrixStruct(testFN);
		MatrixStruct geneLength = new MatrixStruct("");
		test.transpose();
		//test.TPM(geneLengths, true);
		test.transpose();
		test.write(testFN.replace(".txt", "test_TPM.txt"));
	}
	public void print()
	{
		System.out.println("test");
	}
}
