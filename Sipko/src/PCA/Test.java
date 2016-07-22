package PCA;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;
import java.util.stream.Collectors;

import org.apache.commons.math3.stat.StatUtils;

import umcg.genetica.math.stats.Correlation;

public class Test 
{
	
	public static void main(String[] args) throws IOException, InterruptedException
	{
		MatrixStruct matrixStruct = new MatrixStruct("E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/TESTexpression.txt");
		
		double[] s1 = new double[]{ 5,7,9,  8,56,23,1,15.571428571428571,15.571428571428571,15.571428571428571};
		double[] s2 = new double[]{54,1,2,654, 2,32,1,106.57142857142857,106.57142857142857,106.57142857142857};
		
		double avg1 = StatUtils.mean(s1);
		System.out.println("Average1 =" + avg1);
		double avg2 = StatUtils.mean(s2);
		System.out.println("Average2 =" + avg2);
		double test = Math.log(0);
		if(Double.isInfinite(test))
			System.out.println("test =" + test);
		double r = Correlation.correlate(avg1,avg2, s1,s2);
		System.out.println("correlation =" + r);
		//test(matrixStruct);
		
		
		
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
	private static void test(Object matrixStruct) {
		if(matrixStruct instanceof  MatrixStruct)
		{
			System.out.println("yes!");
			MatrixStruct temp = (MatrixStruct)matrixStruct;
			System.out.println(temp.getRowHeaders()[0]);
		}
		
	}
	public void print()
	{
		System.out.println("test");
	}
}
