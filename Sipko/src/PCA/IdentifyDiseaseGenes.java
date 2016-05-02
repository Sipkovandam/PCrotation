package PCA;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import pca.MatrixStruct;

public class IdentifyDiseaseGenes 
{
	static double cutoff = 3;//number of standard deviations that the average expression over the samples has to be away from the
						//average expression of all genes, to be considered a disease gene
	
	public static void main(String[] args) throws IOException
	{
		//Sample in which you wish to find outlier genes
		String sampleFN = "E:/Groningen/Data/PublicSamples/Test13/directPCA_Voom_0.2/18DownSyndrome26Normal2Cancer_counts/PC_1-300_DevidedBySTdevs.txt";
		
		//String sampleFN = "E:/Groningen/Data/Monogenetic_disease_samples_RadBoud/directPCA_Voom_0.2/CountsGENES_Radboud/PC_1-300_DevidedBySTdevs.txt";
		int startCol = 0;//first column to calculate the average from
		int endCol = 9;// last column to calculate the average from
		int[] cols = new int[]{2,13,18};//{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18};
		String writeName = "E:/Groningen/Data/PublicSamples/Test13/directPCA_Voom_0.2/18DownSyndrome26Normal2Cancer_counts/19-04-2016/PC_1-300_OutliersCol0-2_3Stdevs.txt";
		
		if(args.length < 2)
			checkArgs(args);
		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase()){
				case "filename":
					sampleFN = value;
					break;
				case "writename":
					writeName = value;
					break;
				case "diseasesamples":
					endCol = Integer.parseInt(value);
					break;
				case "startcol":
					startCol = Integer.parseInt(value);
					break;
				case "endcol":
					endCol = Integer.parseInt(value);
					break;

				default:
					checkArgs(args);
					System.out.println("Incorrect argument supplied; exiting");
					System.exit(1);
			}
		}
		
		if(writeName == null)
			writeName = sampleFN.replace(".txt", "diseaseGenes_"+endCol+"samples_"+cutoff+"STdevs.txt");
		MatrixStruct samples = new MatrixStruct(sampleFN);
		samples.putGenesOnRows();
		
		//first calculate the average expression and STdev of all genes over the disease samples
		MatrixStruct summaryStats = new MatrixStruct(samples.rows(), 2);
		summaryStats.setRowHeaders(samples.getRowHeaders());
		summaryStats.setColHeaders(new String[]{"Averages","STdevs"});
		
		double[] allAverages = new double[samples.rows()];
		for(int r = 0; r < samples.rows(); r++)
		{
			double[] diseaseValues = new double[cols.length];
			for(int c = 0; c < cols.length; c++)
				diseaseValues[c] = samples.matrix.get(r, cols[c]);
			
			double rowAverage = org.apache.commons.math3.stat.StatUtils.mean(diseaseValues);
			double stdev = Math.pow(org.apache.commons.math3.stat.StatUtils.variance(diseaseValues),0.5);
			
			allAverages[r] = rowAverage;
			summaryStats.matrix.set(r, 0, rowAverage);
			summaryStats.matrix.set(r, 1, stdev);
		}
		//Identify all genes that are 2.5 standard deviations away from the average expression of genes (average will be around 0)
		double average = org.apache.commons.math3.stat.StatUtils.mean(allAverages);
		double stdev = Math.pow(org.apache.commons.math3.stat.StatUtils.populationVariance(allAverages),0.5);
		System.out.println("Average of all genes =" + average + " stdev =" + stdev);
		//ArrayList<String> diseaseGenes = new ArrayList<String>();
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(writeName)));
		writer.write("Disease Genes\tAverage\tstdev\n");
		for(int r = 0; r < summaryStats.rows(); r++)
		{
			double rowAvg = summaryStats.matrix.get(r, 0);
			double rowStdev = summaryStats.matrix.get(r, 1);
			if(rowAvg > average+(cutoff*stdev) || rowAvg< average-(cutoff*stdev))
			{
				//check if the stdev within the disease samples is not to big
				//Identify among that list of genes all genes that have a standard deviation < average
				if(rowStdev < rowAvg || rowStdev < -rowAvg)
				{
					//add it to the list
					//diseaseGenes.add(summaryStats.getRowHeaders()[r]);
					writer.write(summaryStats.getRowHeaders()[r]+"\t" +rowAvg + "\t" + rowStdev + "\n");
				}
			}
				
		}
		writer.close();
		System.out.println("Done! File written to: " + writeName);
	}
	public static void checkArgs(String[] args)
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
			return;
		System.out.println("This script identifies outlier genes across a number of samples\n"
				+ "fileName=<fileName> - FileName containing the corrected samples\n"
				+ "saveName=<savefileName> - Name of the file containing the outlier/disease genes (optional)\n"
				+ "diseasesamples=<number> - Number of disease samples in the file, these have to be the first <number> columns\n");
		System.exit(1);
	}
}
