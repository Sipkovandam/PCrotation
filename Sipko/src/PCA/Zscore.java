package PCA;

import java.io.IOException;
import java.lang.Math;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.stat.inference.WilcoxonSignedRankTest;
import org.apache.commons.math3.stat.inference.TTest;
import org.apache.commons.math3.stat.StatUtils;
import pca.MatrixStruct;


public class Zscore
{
	public static double zScoreCutoff = 2;
	
	public static void main(String[] args) throws IOException
	{
		String sampleScoreFN = "";
		MatrixStruct samples = new MatrixStruct(sampleScoreFN);

		MatrixStruct output = changeToZscores(samples);
		output.write(sampleScoreFN.replace(".txt", "_zScores.txt"));
	}

	public static MatrixStruct changeToZscores(MatrixStruct samples) 
	{
		return changeToZscores(samples, null);
	}
	
	public static MatrixStruct changeToZscores(MatrixStruct samples, String statsFileName) 
	{
		System.out.println(" Memory usage:" + Runtime.getRuntime().totalMemory()/1024/1024/1024 + " GB");
		MatrixStruct scoreStats = new MatrixStruct(samples.rows(),2);
		scoreStats.setRowHeaders(samples.getRowHeaders());
		scoreStats.setColHeaders(new String[]{"Average","StDev"});
		if(statsFileName != null)
			scoreStats = new MatrixStruct(statsFileName);
		
		for(int r = 0; r < samples.getRowHeaders().length; r++)//each has the PCscores for 1 PC for all samples
		{
			double mean = 0,stDev = 0;
			if(statsFileName == null)
			{
				double[] rowValues = samples.getRowValues(r);
				mean = org.apache.commons.math3.stat.StatUtils.mean(rowValues);
				double variance = org.apache.commons.math3.stat.StatUtils.populationVariance(rowValues);
				stDev = java.lang.Math.pow(variance,0.5);
				scoreStats.matrix.set(r, 0, mean);
				scoreStats.matrix.set(r, 1, stDev);
			}
			else
			{
				mean = scoreStats.matrix.get(r, 0);
				stDev = scoreStats.matrix.get(r, 1);
			}
			
			for(int c = 0; c < samples.getColHeaders().length; c++)
			{
				double currentVal = samples.matrix.get(r, c);
				samples.matrix.set(r, c, zScore(currentVal,mean,stDev));
			}
			
//			Arrays.stream(rowValues)
//				  .map(value -> zScore(value,mean,stDev))
//				  .forEach(zScore -> storePCinfo(zScore,output));
			
		}
		return scoreStats;
	}

	private static double zScore(double value, double mean, double stDev) {
		if(stDev == 0) return 0;
		return (value-mean)/stDev;
	}
}

//double zScore = 0;
//new WilcoxonSignedRankTest().wilcoxonSignedRank(arg0, arg1)
//WilcoxonSignedRankTest test = new WilcoxonSignedRankTest();
//TTest tTest = new TTest();
//tTest.
//tTest.
//test.
//samples.getRowValues(r);
