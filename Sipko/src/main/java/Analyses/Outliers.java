package Analyses;

import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Hashtable;

import MatrixScripts.MatrixStruct;
import umcg.genetica.containers.Pair;

public class Outliers {
	//This script identifies genes in sample that have a high z-score in a specific sample but in no others 
	
	//does not identify lower yet, need to add
	
	public static String sampleFN = "E:/Groningen/Data/PublicSamples/09-2016/est_counts_nocancernocellline_Cov_Rlog/clinicSamples/normalized_log2_zScores1.txt";//PC_1-300__zScores.txt";//normalized_log2_zScores1.txt";
	public static String sampleAvgFN = "E:/Groningen/Data/PublicSamples/09-2016/est_counts_nocancernocellline_Cov_Rlog/SAMPLE_Norm_GeneAverages.txt";
	public static String sampleStdevFN = "E:/Groningen/Data/PublicSamples/09-2016/est_counts_nocancernocellline_Cov_Rlog/gene_STDevs.txt";
	public static double minGap = 1;//minimun difference there should be between a z-score of this sample and all others.
	public static double minCOV = 5;//minimum inverse coefficient of variation a gene should have in the public data in order to be considered
	static boolean underexpressed = false; //if true finds the underexpressed outlier, if false the over expressed ones
	
	public static void main (String[] args)
	{
		MatrixStruct samples = new MatrixStruct(sampleFN);
		MatrixStruct avg = new MatrixStruct(sampleAvgFN);
		MatrixStruct stdevs = new MatrixStruct(sampleStdevFN);
		MatrixStruct stats = avg.mergeColumns(stdevs);
		
		samples = getHighCOV(samples, stats);//only keep genes that have a high inverse coefficient of variation in the public data 
		
		samples.putGenesOnRows();
		Hashtable<String, ArrayList<String>> outliersPerSample = getOutliersPerSample(samples);
		
		Enumeration<String> sampleKeys = outliersPerSample.keys();
		while(sampleKeys.hasMoreElements())
		{
			String sample = sampleKeys.nextElement();
			
			String line = sample;
			ArrayList<String> genes = outliersPerSample.get(sample);
			for(int g = 0; g < genes.size(); g++)
				line += "\t"+ genes.get(g);
			//if(sample.contains("160613_SN163_0713_AC8NKTACXX_L5_CACCGG"))
				System.out.println(line);
		}
	}

	private static Hashtable<String, ArrayList<String>> getOutliersPerSample(MatrixStruct samples) {
		Hashtable<String, ArrayList<String>> outliers = new Hashtable<String, ArrayList<String>>();
		for(int r = 0; r< samples.rows(); r++)
		{
			int index = 0;
			Pair<Integer, Integer> pair = null;
			if(!underexpressed)
				pair = getOutlierLargest(samples, r);
			else
				pair = getOutlierSmallest(samples, r);
			
			int n = pair.getLeft();
			index = pair.getRight();
			
			if(n<=1)//if it is a single outlier (thus not also an outlier in any other samples) as defined by being 1 Standard deviation away more then any other sample
			{
				ArrayList<String> sampleOutliers = outliers.get(samples.getColHeaders()[index]);
				if(sampleOutliers == null)
					sampleOutliers = new ArrayList<String>();
				sampleOutliers.add(samples.getRowHeaders()[r]);
				outliers.put(samples.getColHeaders()[index], sampleOutliers);
			}
		}
		return outliers;
	}

	private static Pair<Integer, Integer> getOutlierSmallest(MatrixStruct samples, int r) {
		Pair<Integer,Double> getLarge = samples.getRowSmallest(r);
		double smallest = getLarge.getRight();
		int index = getLarge.getLeft();
		int n = 0;
		for(int c = 0; c < samples.cols(); c++)
		{
			if(smallest > samples.matrix.get(r, c)-minGap)//how often is the smallest value bigger then any other value minus the gap you want to be there
				n++;
		}
		return new Pair<Integer, Integer>(n,index);
	}

	private static Pair<Integer, Integer> getOutlierLargest(MatrixStruct samples, int r) {
		Pair<Integer,Double> getLarge = samples.getRowLargest(r);
		double largest = getLarge.getRight();
		int index = getLarge.getLeft();
		int n = 0;
		for(int c = 0; c < samples.cols(); c++)
		{
			if(largest < samples.matrix.get(r, c)+minGap)
				n++;
		}
		return new Pair<Integer, Integer>(n,index);
	}

	private static MatrixStruct getHighCOV(MatrixStruct samples, MatrixStruct stats) {
		HashMap<String, Integer> keepgenes = new HashMap<String, Integer>();
		int outRow = 0;
		for(int r = 0; r< stats.rows(); r++)
		{
			double cov = stats.matrix.get(r, 0)/stats.matrix.get(r, 1);//=average/stdev
			if(cov > minCOV)
			{
				keepgenes.put(stats.getRowHeaders()[r], outRow);
				outRow++;
			}
		}
		samples.keepIDs(keepgenes);
		
		return samples;
	}
	
}
