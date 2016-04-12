package PCA;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;

import org.apache.commons.math3.stat.inference.TTest;

import pca.MatrixStruct;
import umcg.genetica.math.stats.Correlation;

public class PCcorrection 
{
	
	static boolean setLowestToAverage = false;
	static boolean adjustSampleAverages = true;
	
	public static void main (String[] args) throws IOException
	{
		//String sampleFile = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/" + "RandomSamples.txt";
		String sampleFile = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/" + "TESTexpression.txt";
		String vectorFolder = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/TESTexpression/";
		String writeFolder = null;
		String chr21FN = null;//"E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/Chr21FakeForTest.txt";
		
//		String sampleFile = "E:/Groningen/Data/PublicSamples/Test12/" + "18DownSyndrome26Normal2Cancer_TPM.txt";
//		String vectorFolder = "E:/Groningen/Data/PublicSamples/Test13/TPM_9900SamplesTop0.2/";
//		String writeFolder = "E:/Groningen/Data/PublicSamples/Test13/TPM_9900SamplesTop0.2/18DownSyndrome26Normal2Cancer_TPM_N/";
		
//		String sampleFile = "E:/Groningen/Data/Monogenetic_disease_samples_RadBoud/" + "CountsGENES.txt";
//		String vectorFolder = "E:/Groningen/Data/Monogenetic_disease_samples_RadBoud/directPCA_Rlog_1.0_Correl/";
//		String writeFolder = "E:/Groningen/Data/Monogenetic_disease_samples_RadBoud/directPCA_Rlog_1.0_Correl/RadboudSamples/";
		
//		String sampleFile = "E:/Groningen/Data/Iris/CountsSD200/" + "CountsGENES.txt";
//		String vectorFolder = "E:/Groningen/Data/PublicSamples/04-2016/directPCA_Rlog_0.2/";
//		String writeFolder = "E:/Groningen/Data/Iris/CountsSD200/directPCA_Rlog_0.2_covar/IrisSD200_2/";
		
		String chromLocationsFile = "E:/Groningen/Data/GenePositionInfo_Chr1-22.txt";
		//String chromLocationsFile = "E:/Groningen/Data/GenePositionInfo_23X_24Y_25MT_26rest.txt";
		
		boolean log2 = true;
		boolean correctInputForSTdevs = false;
		boolean correctResultsForSTdevs = true;
		boolean skipQuantileNorm = true;
		
		int optimalPCremoval = -1;
		
		double zScoresCutoff = Double.parseDouble("0");//I removed this (commented out), as it seems pretty useless
		double correctTotalReadCount = -1;//log((gene+0.5)/total*value) //value = number of counts the sample ends up having in total (does noting if 0 or lower)
		double spearman = -1;//-1 is optimal per each sample. Any other "number" than 0 is optimal per "number" samples, where it benefits at least half those "number" of samples
		double addBeforeLog = 0;
		double rLog = 1000000;
		
		String PCs = "1-211";//null if you don't want to set any PCs to correct for (just corrects for the defaults (100,300,500,1000,5000)
		
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
				case "chrom":
					chromLocationsFile = value;
//				case "zscorescutoff":
//					zScoresCutoff = Double.parseDouble(value);
//					break;
				case "correctinputforstdevs":
					correctInputForSTdevs = Boolean.parseBoolean(value);
					break;
				case "log2":
					log2 = Boolean.parseBoolean(value);
					break;
				case "correctresultsforstdevs":
					correctResultsForSTdevs = Boolean.parseBoolean(value);
					break;
				case "noqn":
					skipQuantileNorm = Boolean.parseBoolean(value);
					break;
				case "writefolder":
					writeFolder = value;
					break;
				case "chr21":
					chr21FN = value;
					break;
				case "optimalpcremoval":
					optimalPCremoval = Integer.parseInt(value);
					break;
				case "pcs":
					PCs = value;
					break;
				case "correcttotalreadcount":
					correctTotalReadCount = Double.parseDouble(value);
					break;
				case "spearman":
					spearman = Double.parseDouble(value);
					break;
				case "rlog":
					rLog = Double.parseDouble(value);
					break;
				case "lowesttoaverage":
					setLowestToAverage = Boolean.parseBoolean(value);
					break;
				case "adjustsampleaverages":
					adjustSampleAverages = Boolean.parseBoolean(value);
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
		
		writeParameters(sampleFile, vectorFolder, writeFolder, log2, correctInputForSTdevs, zScoresCutoff, 
				correctResultsForSTdevs, skipQuantileNorm, PCs, correctTotalReadCount, spearman, rLog);
		
		MatrixStruct[] rotationMatrixes = RotateSample.rotate(sampleFile, vectorFolder, writeFolder, correctInputForSTdevs, 
				log2, skipQuantileNorm, correctTotalReadCount, spearman, adjustSampleAverages, setLowestToAverage,
				addBeforeLog, rLog, sampleFile, chromLocationsFile);
		MatrixStruct sampleStruct = rotationMatrixes[2];
		System.out.println("rows = " + sampleStruct.rows());
		
		MatrixStruct zScoreMatrix = null;
//		pca.PCA.log("13. Calculating zScores");
//		String zScoreStats = vectorFolder+"pcZscores_Stats.txt";
//		MatrixStruct zScoreMatrix = rotationMatrixes[0].copy();
//		Zscore.changeToZscores(zScoreMatrix, zScoreStats);
//		
//		pca.PCA.log("14. Writing zScores");
//		zScoreMatrix.write(writeFolder+"pcZscoresSamples.txt");

		String scoreFile = writeFolder + "SAMPLE.PC.scores.txt";
		MatrixStruct scores = new MatrixStruct(scoreFile);
	
		MatrixStruct eigenVectors = rotationMatrixes[3];
		
		pca.PCA.log("16. Adjusting for PCs");
		
		//some very sloppy code. Just want to do things quick atm...
		int[] PCAadjustments = new int[]{0,100,1,25,300,500,1000,5000,eigenVectors.rows()};
		MatrixStruct chr21 = null;
		if(chr21FN != null && new File(chr21FN).exists())
			chr21 = new MatrixStruct(chr21FN);
			
		adjustForPCs(sampleStruct, PCAadjustments, eigenVectors, scores, writeFolder, vectorFolder, 
				zScoreMatrix, zScoresCutoff, log2, correctResultsForSTdevs,chr21, optimalPCremoval, PCs);
		
		System.out.println("Done, Results saved in: " + writeFolder);
	}
	private static void writeParameters(String sampleFile, String vectorFolder, String writeFolder, boolean log2,
			boolean correctInputForSTdevs, double zScoresCutoff, boolean correctResultsForSTdevs, boolean skipQN, String PCs,
			double correctTotalReadCount, double spearman, double rLog) throws IOException {
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
		writer.write("SkipQuantileNormalization\t" + skipQN + "\n");
		writer.write("PCs to correct\t" + PCs + "\n");
		writer.write("correctTotalReadCount\t" + correctTotalReadCount + "\n");
		writer.write("spearman\t" + spearman + "\n");
		writer.write("lowesttoaverage\t" + setLowestToAverage + "\n");
		writer.write("adjustSampleAverages\t" + adjustSampleAverages + "\n");
		writer.write("rLog\t" + rLog + "\n");
		writer.close();
		
	}
	private static void adjustForPCs(MatrixStruct inputMatrix, int[] PCAadjustments, MatrixStruct eigenVectors, MatrixStruct scores, 
			String writeFolder, String vectorFolder, MatrixStruct zScores, double zScoresCutoff, boolean log2, 
			boolean correctResultsForSTdevs, MatrixStruct chr21, int optimalPCremoval, String PCs) throws IOException 
	{
//		pca.PCA.log("Copying matrix");
//		MatrixStruct sampleStruct = inputMatrix.copy();
//		pca.PCA.log("Matrix copy done");
//		if(zScoresCutoff >= 0)//if adjusting based on z-scores get the PCs that should be removed
//		{
//			for(int s = 0; s < sampleStruct.cols();s++)
//			{	
//				ArrayList<Integer> PCsToAdjust = new ArrayList<Integer>();
//				String outputString = "";
//				
//				
//				correctPCs(sampleStruct, scores, eigenVectors, PCsToAdjust, s);
//				
////				if(zScores != null)
////					System.out.println("Sample: "+ s + "/" + zScores.cols()+ "-->" + sampleStruct.getColHeaders()[s] + " Removing " + PCsToAdjust.size() + "/ "+ zScores.rows() +" PCs; Remaining PCs:" + outputString);
//
//			}
//			String writeFileName = writeFolder+"ZscoreCutoff_" + zScoresCutoff + ".txt";
//			sampleStruct.write(writeFileName);
//			createSmoothedFiles(sampleStruct,correctResultsForSTdevs,vectorFolder, writeFileName);
//		}

		MatrixStruct tTestResults = null;
		MatrixStruct difference = null;
		
		ArrayList<Integer> userList = parsePCs(PCs);
		//System.out.println("userList=" + userList);
		pca.PCA.log("Calculating variance explained");
		varianceExplained(PCAadjustments[PCAadjustments.length-1], writeFolder, eigenVectors, inputMatrix);
		pca.PCA.log("Calculating variance explained done");
		
		for(int pcs = 0; pcs < PCAadjustments.length; pcs++)//for all the different numbers of PCs to correct for
		{
			MatrixStruct sampleStruct = inputMatrix.copy();
			int adjustPCs = PCAadjustments[pcs];
			if(adjustPCs > eigenVectors.rows())
				adjustPCs = eigenVectors.rows()+1;
			String writePCName = Integer.toString(adjustPCs);
			
			ArrayList<Integer> PCsToAdjust = new ArrayList<Integer>();
			for(int p = 1; p < adjustPCs; p++) { PCsToAdjust.add(p);}			
			if(pcs == 1 && userList != null)//if the user has a specific list he wants to correct for use this instead of the 100 PC correction.
			{
				PCsToAdjust = userList;
				writePCName = PCs; 
			}
//				if(zScoresCutoff > 0)//for the remaining PCs correct only if Z-score is small
//				{
//					for(int z = adjustPCs; z < zScores.rows(); z++)
//					{
//						double zScore = zScores.matrix.get(z, s);
//						if(zScore<zScoresCutoff && zScore > -zScoresCutoff)//add all PCs < then cuttoff and bigger then -cutoff
//							PCsToAdjust.add(z+1);
//					}
//				}
			if(PCsToAdjust.size() == 0 || chr21 == null)//!= PCAadjustments.length-1
				correctPCs(sampleStruct, scores, eigenVectors, PCsToAdjust,null,null, null, optimalPCremoval);
			else //if it is correcting for the full batch of PCs, also calculate how significant the expression changes are for genes on chr 21 for each PC that is subtracted
			{
				if(chr21 != null)
				{
					chr21.keepRows1Matrix(sampleStruct);
					
					tTestResults = new MatrixStruct(PCsToAdjust.size(), sampleStruct.cols());
					difference = new MatrixStruct(PCsToAdjust.size(), sampleStruct.cols());
					String[] rowHeaders = new String[PCsToAdjust.size()];
					for(int p = 0; p < PCsToAdjust.size(); p++)
					{
						rowHeaders[p] = "PC"+Integer.toString(p+1);
					}
					tTestResults.setRowHeaders(rowHeaders);
					tTestResults.setColHeaders(sampleStruct.getColHeaders());
					difference.setRowHeaders(rowHeaders);
					difference.setColHeaders(sampleStruct.getColHeaders());
					
					if(chr21.rows()>1)
						correctPCs(sampleStruct, scores, eigenVectors, PCsToAdjust,chr21,tTestResults, difference, optimalPCremoval);
					else
						System.out.println("WARNING!: less then 2 genes of chr21 are present and thus no statistical test can be conducted on the difference after the correction of each PC");
				}
			}
			
			String writeFileName = writeFolder+"PC_1-"+writePCName+"_.txt";
			sampleStruct.write(writeFileName);
			createSmoothedFiles(sampleStruct,correctResultsForSTdevs,vectorFolder, writeFileName);
		}
		if(chr21 != null)
		{
			tTestResults.write(writeFolder+"tTestPerPCchr21.txt");
			difference.write(writeFolder+"differencePerPCchr21.txt");
		}
	}
	private static void varianceExplained(int i, String writeFolder, MatrixStruct eigenVectors, MatrixStruct sampleStruct) throws IOException {
		MatrixStruct explained = new MatrixStruct(eigenVectors.rows(), sampleStruct.cols());
		explained.setRowHeaders(eigenVectors.getRowHeaders());
		explained.setColHeaders(sampleStruct.getColHeaders());
		for(int c = 0; c < explained.cols(); c++)//samples are on the columns
			for(int r = 0; r < explained.rows(); r++)//PCs are on the rows
			{
				double[] evValues = eigenVectors.getRowValues(r);
				double[] sampleValues = sampleStruct.getColValues(c);
				double correlation = Correlation.correlate(evValues,sampleValues);
				double varianceExplained = Math.pow(correlation, 2);
				explained.matrix.set(r,c,varianceExplained);
			}
		explained.write(writeFolder+"VarianceExplained.txt");
	}
	private static void correctPCs(MatrixStruct sampleStruct, MatrixStruct scores, MatrixStruct eigenVectors,
			ArrayList<Integer> PCsToAdjust, MatrixStruct chr21, MatrixStruct tTestResults,MatrixStruct difference
			, int optimalPCremoval) {
		double prevPvalue = 1;
		int outcol = 0;
		for(int pc : PCsToAdjust)//correct this sample for the selected PCs
		{
			if(pc > eigenVectors.rows())
				break;
			for(int s = 0; s < sampleStruct.cols();s++)
			{	
				//remove the signal of this single PC from all the genes
				double[][] out = removeSignalAllgenes(sampleStruct, chr21, scores, pc, s, eigenVectors, false);
				double[] onChr21 = out[0];
				double[] others = out[1];
				//check if there is a significant difference
				if(chr21 != null)
				{
					TTest tTest = new TTest();
					double pValue = tTest.tTest(onChr21,others)/2;
					double avgChr21 = org.apache.commons.math3.stat.StatUtils.mean(onChr21);
					double avgOthers = org.apache.commons.math3.stat.StatUtils.mean(others);
					double diff= avgOthers-avgChr21;
					
					if(optimalPCremoval == 1)
						if(prevPvalue < pValue)//if the previous p-value is smaller, just add the signal back on
						{
							removeSignalAllgenes(sampleStruct, chr21, scores, pc, s, eigenVectors, true);//add signal back on
							pValue = prevPvalue;
						}
					prevPvalue = pValue;
					
					tTestResults.matrix.set(outcol, s, pValue);
					difference.matrix.set(outcol, s, diff);
				}
			}
			outcol++;
			if(optimalPCremoval>0 && pc>1)//add signal back on that decreases the signal difference between chr21 and the other genes in at least half the down samples
			{
				//Find out if this last PC should have been kept
				int n = 0;
				for(int c = 0; c < optimalPCremoval; c++)
				{
					if(tTestResults.matrix.get(pc-1-1,c) < tTestResults.matrix.get(pc-1,c))
						n++;
				}
				if(n>= optimalPCremoval/2)//the number of samples in which the PC has to be "bad"(decreasing the difference between genes on chr21 and the other genes)
				{
					for(int s = 0; s < sampleStruct.cols();s++)
					{	
						removeSignalAllgenes(sampleStruct, chr21, scores, pc, s, eigenVectors, true);//add signal back on
					}
				}
			}	
		}	
	}
	private static double[][] removeSignalAllgenes(MatrixStruct sampleStruct, MatrixStruct chr21, MatrixStruct scores, 
			int pc, int s, MatrixStruct eigenVectors, boolean add) {
		//correct all the genes for this PC
		double[] onChr21 = null,others =null;
		//, beforeOthers = null,Others =null
		if(chr21 != null)
		{
			onChr21 = new double[chr21.rows()];
			others = new double[sampleStruct.rows()-chr21.rows()];
		}
		int chr21Index = 0;
		int othersIndex = 0;
		
		for(int gene = 0; gene < sampleStruct.rows(); gene++)
		{
			double signal = scores.matrix.get(pc-1, s)*eigenVectors.matrix.get(pc-1, gene);
//			if(gene == 0 && pc <= 2)
//				System.out.println("gene = " +sampleStruct.getRowHeaders()[gene] + " PC =" +pc+" signal = \t" + signal + 
//						" score ="+ scores.matrix.get(pc-1, s) + " vectorValue= " + eigenVectors.matrix.get(pc-1, gene)+ "val before = " +sampleStruct.matrix.get(gene, s));
			if(add)//add the signal instead of removing it
				signal *= -1;
			sampleStruct.matrix.add(gene, s, -signal);
//			if(gene == 0 && pc <= 2)
//				System.out.println("gene = " +sampleStruct.getRowHeaders()[gene] + " PC =" +pc+"val after = " +sampleStruct.matrix.get(gene, s));
			if(chr21 != null)
			{
				if(chr21.rowHash.containsKey(sampleStruct.getRowHeaders()[gene]))
				{
					onChr21[chr21Index] = sampleStruct.matrix.get(gene,s);
					chr21Index++;
				}
				else
				{
					others[othersIndex] = sampleStruct.matrix.get(gene,s);
					othersIndex++;
				}
			}
		}
		return new double[][]{onChr21,others};
	}
	private static void createSmoothedFiles(MatrixStruct sampleStruct, boolean correctResultsForSTdevs,
			String vectorFolder, String writeFileName) throws IOException {

		devideBySTdev(sampleStruct,correctResultsForSTdevs,vectorFolder, writeFileName);
		
		int smoothNgenes = 10;
//		pca.PCA.log("18. Smooth "+smoothNgenes+" genes");
//		smoothSignal(sampleStruct.copy(),smoothNgenes, writeFileName);
		
		smoothNgenes = 100;
		pca.PCA.log("19. Smooth "+smoothNgenes+" genes");
		smoothSignal(sampleStruct.copy(),smoothNgenes,writeFileName);
		
//		double topPercent = 0.5;
//		keepTopPercentage(sampleStruct, vectorFolder+"SAMPLE_Norm_columnAverages.txt",topPercent, writeFileName, false);
//		
//		smoothNgenes = 10;
//		pca.PCA.log("21. Smooth "+smoothNgenes+" genes");
//		smoothSignal(sampleStruct.copy(),smoothNgenes, writeFileName.replace(".txt", "Top"+topPercent*100+"%_.txt"));
//		
//		smoothNgenes = 100;
//		pca.PCA.log("22. Smooth "+smoothNgenes+" genes");
//		smoothSignal(sampleStruct.copy(),smoothNgenes,writeFileName.replace(".txt", "Top"+topPercent*100+"%_.txt"));
		
	}
	static void keepTopPercentage(MatrixStruct sampleStruct, String averagesFN, double topPercent, String writeFileName, boolean lowest, boolean writeAll) throws IOException {
		pca.PCA.log("20. Highest "+ topPercent*100 + "% only");
		MatrixStruct averages = new MatrixStruct(averagesFN);
		averages.sortCol(0);
		averages.write(averagesFN.replace(".txt", "SAMPLE_Norm_columnAverages_SORTED_ALL.txt"));
		MatrixStruct part = keepPart(averages, topPercent, lowest);
		averages.write(averagesFN.replace(".txt", "SAMPLE_Norm_columnAverages_SORTED.txt"));
		sampleStruct.keepRows(averages);
		part.write(averagesFN.replace(".txt", "SAMPLE_Norm_columnAverages_"+topPercent+"highestExpressed.txt"));
		sampleStruct.keepRows(part);
		if(writeAll)
			sampleStruct.write(writeFileName.replace(".txt",".Top"+topPercent*100+"%only.txt"));
	}
	private static void devideBySTdev(MatrixStruct sampleStruct, boolean correctResultsForSTdevs, String vectorFolder,
			String writeFileName) throws IOException {

		if(correctResultsForSTdevs)
		{
			pca.PCA.log("16. Divide by standard deviation");
			MatrixStruct STdevs = new MatrixStruct(vectorFolder+"gene_STDevs.txt");
			//if the standard deviation is smaller then 1, set it to 1 to avoid inflated values for genes that have a very small stdev
			for(int s = 0; s < STdevs.rows(); s++)
				if(STdevs.matrix.get(s, 0) < 1)
					STdevs.matrix.set(s, 0,1);
			sampleStruct.divideBy(STdevs, true);
			sampleStruct.write(writeFileName.replace(".txt", "DevidedBySTdevs.txt"));
		}
		
	}
	private static MatrixStruct keepPart(MatrixStruct averages, double topPercent, boolean lowest) //returns the last part of the matrix
	{
		int topX = (int)(topPercent*averages.rows());
		MatrixStruct part = new MatrixStruct(topX, averages.cols());
		part.setColHeaders(averages.getColHeaders());
		if(!lowest)
			for(int r = 0; r < topX; r++)
			{
				part.setRow(r, averages.getRowHeaders()[r], averages.getRowValues(r));
			}
		else
			for(int r = averages.rows()-topX, out = 0; r < averages.rows(); r++, out++)
			{
				part.setRow(out, averages.getRowHeaders()[r], averages.getRowValues(r));
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
	public static ArrayList<Integer> parsePCs(String PCsToAdjust) 
	{//function takes format like: "1,4,5-10,3-10"
		if(PCsToAdjust.contains("null"))
			return null;
		ArrayList<Integer> PCs = new ArrayList<Integer>();
		if(PCsToAdjust.contains("-"))
		{
			String[] eles = PCsToAdjust.split("-");
			for(int e = 0; e < eles.length-1; e++)
			{
				String[] ele = eles[e].split(",");
				int last = ele.length-1;
				int start = Integer.parseInt(ele[last]);
				String[] ele2 = eles[e+1].split(",");
				int end = Integer.parseInt(ele2[0]);
				for(int n = start; n <= end; n++)
					PCs.add(n);
			}
		}
		if(PCsToAdjust.contains(","))
		{
			String[] eles = PCsToAdjust.split(",");
			for(int e = 0; e < eles.length; e++)
			{
				if(!eles[e].contains("-"))
					PCs.add(Integer.parseInt(eles[e]));
			}
		}
		
//		System.out.println("size = " + PCs.size());
//		for(int PC : PCs)
//		{
//			System.out.println("pc = " + PC);
//		}
		
		return PCs;
	}
	public static void checkArgs(String[] args)
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
			return;
		System.out.println("This script calculates the eigenvectors over the genes and uses the following input:\n"
				+ "fileName=<fileName> - Expression file (samples on rows, gene names on columns)\n"
				+ "vectorFolder=<vectorfolder> - Folder with the gene eigen vectors created by CreageGeneEigenvectorFile.java.\n"
				+ "correctinputforstdevs=<false/true> - Corrects for STdevs prior to centering the input data (default=false)\n"
				+ "log2=<false/true> - Log2 transformation after quantile normalization (default=true)\n"
				+ "correctresultsforstdevs=<false/true>  - Divide values obtained after correcting the data by standard deviation (default = true)\n"
				+ "tpm=<false/true>  - Input values are TPM, make sure the eigenvector matrix was also calculated based on TPM values (default = false)\n"
				+ "chr21FN=<fileName>  - File with geneNames. Script will calculate the p-value between the difference between genes"
				+ "\tin this file and all remaining genes (default = null)\n"
				+ "writeFolder=<folderName> - folder to write the results in (default=filename.replace(.txt,/)\n"
				+ "PCs=<PCs> - PCs to correct for. The following formats can be used: "
					+ "   1-100"
					+ "   1,2,6,8"
					+ "   or a combination: 1,5-10,66,100-200");
		System.exit(1);
	}
}

//case "filename":
//	sampleFile = value;
//	break;
//case "vectorfolder":
//	vectorFolder = value;
//	break;

//case "correctinputforstdevs":
//	correctInputForSTdevs = Boolean.parseBoolean(value);
//	break;
//case "log2":
//	log2 = Boolean.parseBoolean(value);
//	break;
//case "correctresultsforstdevs":
//	correctResultsForSTdevs = Boolean.parseBoolean(value);
//	break;
//case "tpm":
//	tpm = Boolean.parseBoolean(value);
//	break;
//case "writefolder":
//	writeFolder = value;
//	break;
//case "chr21":
//	chr21FN = value;
//	break;
//case "optimalpcremoval":
//	optimalPCremoval = Integer.parseInt(value);
//	break;