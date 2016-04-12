package PCA;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;
import java.util.stream.Collectors;

import pca.MatrixStruct;

public class Test 
{
	
	public static void main(String[] args) throws IOException, InterruptedException
	{
		String nan = "NaN";
		double parsed = Double.parseDouble(nan);
		System.out.println(parsed);
		if(Double.isNaN(parsed))
			System.out.println("yes");
		
		/*
		ArrayList<Double> vals = new ArrayList<Double>();
		vals.add(1.0);
		vals.add(11.0);
		vals.add(8.0);
		vals.add(1.0);
		vals.add(111.0);
		Collections.sort(vals);
		for(int v = 0; v < vals.size(); v++)
			System.out.println(vals.get(v));
		
		double test = Double.NaN;
		
		if(Double.isNaN(test))
			System.out.println("Hijs is gek he!");
		
		double[] values = new double[]{90.69690754,71.40272395,87.21937366};
		double standardDev = java.lang.Math.pow(org.apache.commons.math3.stat.StatUtils.populationVariance(values),0.5);
		System.out.println(standardDev);

		Random rnd = new Random();
		for(int x = 0; x < 1000; x++)
		{
			System.out.println(rnd.nextGaussian());
		}*/
		//String FN = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/TESTexpression.txt";
		//MatrixStruct test = new MatrixStruct(FN,3,-1);
		
		//test.write(FN.replace(".txt", "_subsection.txt"));
	}
	public void print()
	{
		System.out.println("test");
	}
}
