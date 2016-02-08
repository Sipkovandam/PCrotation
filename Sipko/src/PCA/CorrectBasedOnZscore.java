package PCA;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;

import pca.MatrixStruct;

public class CorrectBasedOnZscore 
{

	public static void main (String[] args) throws IOException
	{
		String sampleFile = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/" + "RandomSamples.txt";
		String vectorFolder = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/TESTexpression/";
		
//		String sampleFile = "E:/Groningen/Data/PublicSamples/Test8/" + "4DownSyndrome3Normal3Cancer_countsNoDupsWithlog.txt";
//		String vectorFolder = "E:/Groningen/Data/PublicSamples/Test8/NODUPS_WITHLOG_NOSTDEV/";
		
		double zScoresCutoff = Double.parseDouble("1");
		String writeFolder = null;
		
		boolean log2 = true;
		boolean correctInputForSTdevs = false;
		boolean correctResultsForSTdevs = true;
		
		if(args.length==0) 
			checkArgs(args);

		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase()){
				case "filename":
					sampleFile = value;
					break;
				case "vectorfolder":
					vectorFolder = value;
					break;
				case "zscorescutoff":
					zScoresCutoff = Double.parseDouble(value);
					break;
				case "correctinputforstdevs":
					correctInputForSTdevs = Boolean.parseBoolean(value);
					break;
				case "log2":
					log2 = Boolean.parseBoolean(value);
					break;
				case "correctresultsforstdevs":
					correctResultsForSTdevs = Boolean.parseBoolean(value);
					break;
				case "writefolder":
					writeFolder = value;
					break;
//				case "duplicate":
//					duplicateCutoff = Double.parseDouble(value);
//					break;
				default:
					checkArgs(args);
					System.out.println("Incorrect argument supplied; exiting");
					System.exit(1);
			}
		}
		if(writeFolder == null)
			writeFolder = sampleFile.replace(".txt", "_Adj/");
		
		writeParameters(sampleFile, vectorFolder, writeFolder, log2, correctInputForSTdevs, zScoresCutoff, correctResultsForSTdevs);
		
		MatrixStruct[] rotationMatrixes = RotateSample.rotate(sampleFile, vectorFolder, writeFolder, correctInputForSTdevs, log2);//1 is centered matrix
		MatrixStruct sampleStruct = rotationMatrixes[2];
		
		pca.PCA.log("13. Calculating zScores");
		String zScoreStats = vectorFolder+"pcZscores_Stats.txt";
		MatrixStruct zScoreMatrix = rotationMatrixes[0].copy();
		Zscore.changeToZscores(zScoreMatrix, zScoreStats);
		
		pca.PCA.log("14. Writing zScores");
		zScoreMatrix.write(writeFolder+"pcZscoresSamples.txt");

		String scoreFile = writeFolder + "SAMPLE.PC.scores.txt";
		MatrixStruct scores = new MatrixStruct(scoreFile);
		
		pca.PCA.log("15. Reading eigenvector matrix");
		MatrixStruct eigenVectors = rotationMatrixes[3];
		
		pca.PCA.log("16. Adjusting for PCs");
		
		//some very sloppy code. Just want to do things quick atm...
		int[] PCAadjustments = new int[]{0,100,265,300,1000,5000};
		adjustForPCs(sampleStruct.copy(), PCAadjustments, eigenVectors, scores, writeFolder, vectorFolder, zScoreMatrix, zScoresCutoff, log2, correctResultsForSTdevs);
		
		System.out.println("Done");
	}
	private static void writeParameters(String sampleFile, String vectorFolder, String writeFolder, boolean log2,
			boolean correctInputForSTdevs, double zScoresCutoff, boolean correctResultsForSTdevs) throws IOException {
		CreateGeneEigenvectorFile.makeFolder(writeFolder);
		DateFormat dateFormat = new SimpleDateFormat("dd/MM/yyyy HH:mm:ss");
		Date date = new Date();
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(writeFolder+"parameters.txt")));
		writer.write("Date\t" + dateFormat.format(date) + "\n");
		writer.write("Input file\t" + sampleFile + "\n");
		writer.write("writeFolder\t" + writeFolder + "\n");
		writer.write("vectorFolder\t" + vectorFolder + "\n");
		writer.write("log2\t" + log2 + "\n");
		writer.write("correctInputForSTdevs\t" + correctInputForSTdevs + "\n");
		writer.write("zScoresCutoff\t" + zScoresCutoff + "\n");
		writer.write("correctResultsForSTdevs\t" + correctResultsForSTdevs + "\n");
		writer.close();
		
	}
	private static void adjustForPCs(MatrixStruct inputMatrix, int[] PCAadjustments, MatrixStruct eigenVectors, MatrixStruct scores, String writeFolder, String vectorFolder, MatrixStruct zScores, double zScoresCutoff, boolean log2, boolean correctResultsForSTdevs) throws IOException 
	{
		MatrixStruct sampleStruct = inputMatrix.copy();
		if(zScoresCutoff >= 0)//if adjusting based on z-scores get the PCs that should be removed
		{
			for(int s = 0; s < sampleStruct.cols();s++)
			{	
				ArrayList<Integer> PCsToAdjust = new ArrayList<Integer>();
				String outputString = "";
				
				for(int z = 0; z < zScores.rows(); z++)
				{
					double zScore = zScores.matrix.get(z, s);
					if(zScore<zScoresCutoff && zScore > -zScoresCutoff)//add all PCs < then cuttoff and bigger then -cutoff
						PCsToAdjust.add(z+1);
					else
						outputString += "_"+(z+1);
				}
				correctPCs(sampleStruct, scores, eigenVectors, PCsToAdjust, s);
				
				if(zScores != null)
					System.out.println("Sample: "+ s + "/" + zScores.cols()+ "-->" + sampleStruct.getColHeaders()[s] + " Removing " + PCsToAdjust.size() + "/ "+ zScores.rows() +" PCs; Remaining PCs:" + outputString);

			}
			String writeFileName = writeFolder+"ZscoreCutoff_" + zScoresCutoff + ".txt";
			sampleStruct.write(writeFileName);
			createSmoothedFiles(sampleStruct,correctResultsForSTdevs,vectorFolder, writeFileName);
		}
		for(int pcs = 0; pcs < PCAadjustments.length; pcs++)//for all the different numbers of PCs to correct for
		{
			sampleStruct = inputMatrix.copy();
			int adjustPCs = PCAadjustments[pcs];
			System.out.println("AdjustPCs =" + adjustPCs);
			if(adjustPCs > eigenVectors.rows())
				return;
			for(int s = 0; s < sampleStruct.cols();s++)
			{		
				ArrayList<Integer> PCsToAdjust = new ArrayList<Integer>();
				PCsToAdjust = new ArrayList<Integer>();
				for(int p = 1; p < adjustPCs; p++) { PCsToAdjust.add(p);}			
				
				correctPCs(sampleStruct, scores, eigenVectors, PCsToAdjust, s);
			}
			String writeFileName = writeFolder+"PC_1-"+adjustPCs+"_.txt";
			sampleStruct.write(writeFileName);
			createSmoothedFiles(sampleStruct,correctResultsForSTdevs,vectorFolder, writeFileName);
		}
	}
	private static void correctPCs(MatrixStruct sampleStruct, MatrixStruct scores, MatrixStruct eigenVectors,
			ArrayList<Integer> PCsToAdjust, int s) {
		for(int pc : PCsToAdjust)//correct this sample for the selected PCs
		{
			//correct all the genes for this PC
			for(int gene = 0; gene < sampleStruct.rows(); gene++)
			{
				double signal = scores.matrix.get(pc-1, s)*eigenVectors.matrix.get(pc-1, gene);
				sampleStruct.matrix.add(gene, s, -signal);
			}
			
		}	
	}
	private static void createSmoothedFiles(MatrixStruct sampleStruct, boolean correctResultsForSTdevs,
			String vectorFolder, String writeFileName) throws IOException {

		devideBySTdev(sampleStruct,correctResultsForSTdevs,vectorFolder, writeFileName);
		
		int smoothNgenes = 10;
		pca.PCA.log("18. Smooth "+smoothNgenes+" genes");
		smoothSignal(sampleStruct.copy(),smoothNgenes, writeFileName);
		
		smoothNgenes = 100;
		pca.PCA.log("19. Smooth "+smoothNgenes+" genes");
		smoothSignal(sampleStruct.copy(),smoothNgenes,writeFileName);
		
		double topPercent = 0.5;
		keepTopPercentage(sampleStruct, vectorFolder+"SAMPLE_QuantNorm_columnAverages.txt",topPercent, writeFileName);
		
		smoothNgenes = 10;
		pca.PCA.log("21. Smooth "+smoothNgenes+" genes");
		smoothSignal(sampleStruct.copy(),smoothNgenes, writeFileName.replace(".txt", "Top"+topPercent*100+"%_.txt"));
		
		smoothNgenes = 100;
		pca.PCA.log("22. Smooth "+smoothNgenes+" genes");
		smoothSignal(sampleStruct.copy(),smoothNgenes,writeFileName.replace(".txt", "Top"+topPercent*100+"%_.txt"));
		
	}
	static void keepTopPercentage(MatrixStruct sampleStruct, String averagesFN, double topPercent, String writeFileName) throws IOException {
		pca.PCA.log("20. Highest "+ topPercent*100 + "% only");
		MatrixStruct averages = new MatrixStruct(averagesFN);
		averages.sortCol(0);
		MatrixStruct part = keepPart(averages, topPercent);
		averages.write(averagesFN.replace(".txt", "SAMPLE_QuantNorm_columnAverages_SORTED.txt"));
		sampleStruct.keepRows(averages);
		part.write(averagesFN.replace(".txt", "SAMPLE_QuantNorm_columnAverages_"+topPercent+"highestExpressed.txt"));
		sampleStruct.keepRows(part);
		sampleStruct.write(writeFileName.replace(".txt",".Top"+topPercent*100+"%only.txt"));
	}
	private static void devideBySTdev(MatrixStruct sampleStruct, boolean correctResultsForSTdevs, String vectorFolder,
			String writeFileName) throws IOException {

		if(correctResultsForSTdevs)
		{
			pca.PCA.log("16. Divide by standard deviation");
			MatrixStruct STdevs = new MatrixStruct(vectorFolder+"gene_STDevs.txt");
			sampleStruct.divideBy(STdevs, true);
			sampleStruct.write(writeFileName.replace(".txt", "DevidedBySTdevs.txt"));
		}
		
	}
	private static MatrixStruct keepPart(MatrixStruct averages, double topPercent) //returns the last part of the matrix
	{
		int topX = (int)(topPercent*averages.rows());
		MatrixStruct part = new MatrixStruct(topX, averages.cols());
		part.setColHeaders(averages.getColHeaders());
		int out = 0;
		for(int r = 0; r < topX; r++)
		{
			part.setRow(r, averages.getRowHeaders()[r], averages.getRowValues(r));
		}
		part.rowHash = MatrixStruct.makeHash(part.getRowHeaders());
		return part;
	}
	private static void smoothSignal(MatrixStruct sampleStruct, int n, String writeFileName) throws IOException 
	{
		if(n > sampleStruct.rows())
			return;
		if(n%2 ==0)//want an uneven number so you can take half before half after + the one itself
		{
			n++;
			//System.out.println("Taking average of " + n + " numbers");
		}
		for(int c = 0; c < sampleStruct.cols(); c++)
		{
			double runningSum = 0;
			double runningAvg = 0;
			double[] avgs = new double[sampleStruct.rows()];
			int count = 0;
		
			for(int r = 0; r< sampleStruct.rows(); r++)
			{
				if(count < n)//for the first n numbers
				{
					runningSum+= sampleStruct.matrix.get(r,c);
					count++;
				}
				else
				{
					runningSum-= sampleStruct.matrix.get(r-n, c);
					runningSum+= sampleStruct.matrix.get(r,c);
				}
				runningAvg= runningSum/count;
				if(r-n/2>=0)
					avgs[r-n/2] = runningAvg;
			}
			//last n/2 rows
			for(int r = sampleStruct.rows(); r< sampleStruct.rows()+n/2; r++)
			{
				runningSum-=sampleStruct.matrix.get(r-n, c);
				count--;
				runningAvg= runningSum/count;
				avgs[r-n/2] = runningAvg;
			}
			sampleStruct.setCol(avgs, c);
		}
		sampleStruct.write(writeFileName.replace(".txt", ".Smoothed"+n+"Genes.txt"));
	}
	public static void checkArgs(String[] args)
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
			return;
		System.out.println("This script calculates the eigenvectors over the genes and uses the following input:\n"
				+ "fileName=<fileName> - Expression file (samples on rows, gene names on columns)\n"
				+ "vectorFolder=<vectorfolder> - Folder with the gene eigen vectors created by CreageGeneEigenvectorFile.java.\n"
				+ "zScoresCutoff=<number> - all PCs under this treshold will be removed (default=-1 (indicating this is not used))\n"
				+ "correctinputforstdevs=<false/true> - Corrects for STdevs prior to centering the input data (default=false)\n"
				+ "log2=<false/true> - Log2 transformation after quantile normalization (default=true)\n"
				+ "correctresultsforstdevs=<false/true>  - Devide values obtained after correcting the data by standard deviation (default = true)\n"
				+ "writeFolder=<folderName> - folder to write the results in (default=filename.replace(.txt,/)");
		System.exit(1);
	}
}