package PCA;

import java.io.IOException;
import java.util.Arrays;
import java.util.Random;
import java.util.stream.Collectors;

import pca.MatrixStruct;

public class Test 
{
	
	public static void main(String[] args) throws IOException, InterruptedException
	{
		double[] values = new double[]{90.69690754,71.40272395,87.21937366};
		double standardDev = java.lang.Math.pow(org.apache.commons.math3.stat.StatUtils.populationVariance(values),0.5);
		System.out.println(standardDev);

		Random rnd = new Random();
		for(int x = 0; x < 1000; x++)
		{
			System.out.println(rnd.nextGaussian());
		}
		//String FN = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/TESTexpression.txt";
		//MatrixStruct test = new MatrixStruct(FN,3,-1);
		
		//test.write(FN.replace(".txt", "_subsection.txt"));
	}
	public void print()
	{
		System.out.println("test");
	}
}
