package PCA;

import java.io.IOException;
import java.lang.Math;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.stat.inference.WilcoxonSignedRankTest;

import Tools.FileUtils;

import org.apache.commons.math3.stat.inference.TTest;
import org.apache.commons.math3.stat.StatUtils;


public class Zscore
{
	//Converts values per row (gene) to Z-scores per gene across all samples
	
	//public static double zScoreCutoff = 2;
	//public static String sampleFN = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/TESTexpression.txt";
	public static String sampleFN = "E:/Groningen/Data/PublicSamples/09-2016/est_counts_nocancernocellline_Cov_Rlog/clinicSamples/PC_1-300_.txt";
	public static String writeFN = null;
	public static String statsFN = "E:/Groningen/Data/PublicSamples/09-2016/est_counts_nocancernocellline_Cov_Rlog/est_counts_nocancernocellline/PC_1-300__Gene_stats.txt";//"E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/TESTexpression_stats.txt";
	static boolean genes = true;
	
	public static void main(String[] args) throws IOException
	{
		checkArgs(args);
		MatrixStruct samples = new MatrixStruct(sampleFN);
		if(genes)
			samples.putGenesOnRows();
		MatrixStruct output = changeToZscores(samples, statsFN);
		if(writeFN==null)
		{
			writeFN= FileUtils.replaceEnd(sampleFN, "_zScores.txt");
		}
		String statsFN= FileUtils.replaceEnd(sampleFN, "_stats.txt");
		samples.write(writeFN);
		output.write(statsFN);
		System.out.println("Done! File written to: " + writeFN);
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
		
		for(int r = 0; r < samples.getRowHeaders().length; r++)//each has the values for 1 row (for all samples)
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
				String geneName = samples.getRowHeaders()[r];
				if(!scoreStats.getRowHash().containsKey(geneName))
					System.out.println("Warning this gene is missing from the stats file:" + geneName);
				int scoreRow = scoreStats.getRowHash().get(geneName);
				mean = scoreStats.matrix.get(scoreRow, 0);
				stDev = scoreStats.matrix.get(scoreRow, 1);
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
	
	static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko") && args.length < 1)
			return;
		if(args.length < 1)
		{
			System.out.println("Script requires the following argumetns:\n"
					+ "1. samplefn=<countfilename.txt> - File with the counts for every sample\n"
					+ "2. genes=<transcripttogenefile.txt> - If true will automatically put genes on rows to calculate z-scores per gene\n"
					+ "3. statsFN=<statsfilename.txt> - File with averages that should be used in column 1 and stdevs in column 2\n"
					+ "4. writeFN=<writename.txt> - Name of the file to be written (default=<sampleFN>+_zScores.txt");
			System.exit(1);
		}
		
		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase())
			{
//				var = new JSONutil<Vars>().read(var.JSON_FN, var);
				case "samplefn":
					sampleFN =value;
					break;
				case "genes":
					genes = Boolean.parseBoolean(value);
					break;
				case "writefn":
					writeFN=value;
					break;
				case "statsfn":
					statsFN=value;
					break;
				default:
					System.out.println("Incorrect argument supplied:\n"+ args[a] +"\nexiting");
					System.exit(1);
			}
		}
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
