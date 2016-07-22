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
	static Var var = new Var();
	
	public static void main (String[] args) throws IOException
	{
		checkArgs(args);
		
		if(var.writeFolderCorrected == null || var.writeFolderCorrected.length() == 0)
		{
			File sample = new File(var.sampleFile);
			var.writeFolderCorrected = var.getFolderName(var.writeFolder)+sample.getName().replace(".txt", "").replace(".gz", "")+"/";
		}
		
		makeFolder(var.writeFolderCorrected);
		var.writeVars(var.writeFolderCorrected+"variables.txt");
		
		//put the matrix in the same space and calculate the PC scores for the PCs defined based on the pulbic data
		MatrixStruct[] rotationMatrixes = RotateSample.rotate(var.sampleFile, var.writeFolder, var.writeFolderCorrected, var.correctInputForSTdevs, 
				var.log2, var.skipQuantileNorm, var.correctTotalReadCount, var.spearman, var.adjustSampleAverages, var.setLowestToAverage,
				var.addLogVal, var.rLog, var.sampleFile, var.chromLocationsFile, var.correctGCSamples);
		MatrixStruct sampleStruct = rotationMatrixes[2];
		System.out.println("rows = " + sampleStruct.rows());
		
		MatrixStruct zScoreMatrix = null;
//		pca.PCA.log("13. Calculating zScores");
//		String zScoreStats = vectorFolder+"pcZscores_Stats.txt";
//		MatrixStruct zScoreMatrix = rotationMatrixes[0].copy();
//		Zscore.changeToZscores(zScoreMatrix, zScoreStats);
//		
//		pca.PCA.log("14. Writing zScores");
//		zScoreMatrix.write(writeFolder+"pcZscoresSamples.txt");\\\\\

		String scoreFile = var.writeFolderCorrected + "SAMPLE.PC.scores.txt";
		MatrixStruct scores = new MatrixStruct(scoreFile);
	
		MatrixStruct eigenVectors = rotationMatrixes[3];
		
		pca.PCA.log("16. Adjusting for PCs");
		
		int[] PCAadjustments = new int[]{0,100,2,25,300,500,1000};//,5000,eigenVectors.rows()};
		MatrixStruct chr21 = null;
		if(var.chr21FN != null && new File(var.chr21FN).exists())
			chr21 = new MatrixStruct(var.chr21FN);
		
		adjustForPCs(sampleStruct, PCAadjustments, eigenVectors, scores, var.writeFolderCorrected, var.writeFolder, 
				zScoreMatrix, var.zScoresCutoff, var.log2, var.correctResultsForSTdevs, chr21, var.optimalPCremoval, var.PCs);
		
		System.out.println("Done, Results saved in: " + var.writeFolderCorrected);
	}
	private static void adjustForPCs(MatrixStruct inputMatrix, int[] PCAadjustments, MatrixStruct eigenVectors, MatrixStruct scores, 
			String writeFolder, String vectorFolder, MatrixStruct zScores, double zScoresCutoff, boolean log2, 
			boolean correctResultsForSTdevs, MatrixStruct chr21, int optimalPCremoval, String PCs) throws IOException 
	{
		//adjustOnZscores();
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

			if(PCsToAdjust.size() == 0 || chr21 == null)//!= PCAadjustments.length-1
				correctPCs(sampleStruct, scores, eigenVectors, PCsToAdjust,null,null, null, optimalPCremoval);
			else //if it is correcting for the full batch of PCs, also calculate how significant the expression changes are for genes on chr 21 for each PC that is subtracted
				correctPCsTrackChr21(sampleStruct,scores, eigenVectors,PCsToAdjust,chr21, tTestResults, difference, optimalPCremoval);
			
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
	private static void correctPCsTrackChr21(MatrixStruct sampleStruct, MatrixStruct scores, MatrixStruct eigenVectors, 
			ArrayList<Integer> PCsToAdjust, MatrixStruct chr21, MatrixStruct tTestResults, MatrixStruct difference, int optimalPCremoval) {
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
	private static void adjustOnZscores() {
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
				
				//check if there is a significant difference between chr21 and the rest
				trackChr21(out, chr21, sampleStruct, scores, pc, s, eigenVectors, optimalPCremoval, prevPvalue, difference, difference, tTestResults, outcol);
			}
			outcol++;
			
			if(optimalPCremoval>0 && pc>1)//add signal back on that decreases the signal difference between chr21 and the other genes in at least half the down samples
				restoreGoodSignal(sampleStruct, chr21, scores, pc, eigenVectors, true, optimalPCremoval, tTestResults);
		}	
	}
	private static void restoreGoodSignal(MatrixStruct sampleStruct, MatrixStruct chr21, MatrixStruct scores, int pc,
			MatrixStruct eigenVectors, boolean b, int optimalPCremoval, MatrixStruct tTestResults) {
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
	
	private static void trackChr21(double[][] out, MatrixStruct chr21, MatrixStruct sampleStruct, MatrixStruct scores, int pc, int s,
			MatrixStruct eigenVectors, int optimalPCremoval, double prevPvalue, MatrixStruct difference, MatrixStruct difference2, MatrixStruct tTestResults, int outcol) 
	{
		double[] onChr21 = out[0];
		double[] others = out[1];
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
			//if(gene == 0 && pc <= 2 && s==0)
//			if(sampleStruct.getRowHeaders()[gene].contains("ENSG00000250360") && pc <= 2 && s==0)
//				System.out.println("gene = " +sampleStruct.getRowHeaders()[gene] + " geneEigen = " + eigenVectors.getColHeaders()[gene] + " PC =" +pc+" signal = \t" + signal + 
//						" score ="+ scores.matrix.get(pc-1, s) + " vectorValue= " + eigenVectors.matrix.get(pc-1, gene)+ "val before = " +sampleStruct.matrix.get(gene, s));
			if(add)//add the signal instead of removing it
				signal *= -1;
			sampleStruct.matrix.add(gene, s, -signal);
			//if(gene == 0 && pc <= 2 && s==0)
//			if(sampleStruct.getRowHeaders()[gene].contains("ENSG00000250360") && pc <= 2 && s==0)
//				System.out.println("gene = " +sampleStruct.getRowHeaders()[gene] + " geneEigen = " + eigenVectors.getColHeaders()[gene] + " PC =" +pc+"val after = " +sampleStruct.matrix.get(gene, s));
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
			STdevs.keepRows(sampleStruct);
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
	static void makeFolder(String writeFolder) 
	{
		File folder = new File(writeFolder);
		if(!folder.exists())
		{
			folder.mkdir();
		}
		
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
	static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
			return;
		if(args.length < 1)
		{
			System.out.println("Script requires the following arguments:\n"
					+ "1. filename=<expressionFN.txt> - Expression file with genes on rows samples on columns\n"
					+ "2. writeFolder=<writeFolderFN.txt> - Folder where the files will be written (default=parentFolder(input.txt))\n"
					+ "3. geoFN=<geoFn.txt> - Optional file with geometric mean per gene to use (default=null)\n");
			System.exit(1);
		}
		
		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase())
			{
				case "json":
					var = var.readVars(value);
				break;
				case "correctresultsforstdevs":
					var.correctResultsForSTdevs = Boolean.parseBoolean(value);
					break;
				case "chr21":
					var.chr21FN = value;
					break;
				case "optimalpcremoval":
					var.optimalPCremoval = Integer.parseInt(value);
					break;
				case "pcs":
					var.PCs = value;
					break;
				case "filename":
					var.sampleFile = value;
					break;
				case "vectorfolder":
					var.writeFolder = value;
					break;
				case "chrom":
					var.chromLocationsFile = value;
				case "correctinputforstdevs":
					var.correctInputForSTdevs = Boolean.parseBoolean(value);
					break;
				case "log2":
					var.log2 = Boolean.parseBoolean(value);
					break;		
				case "noqn":
					var.skipQuantileNorm = Boolean.parseBoolean(value);
					break;
				case "writefolder":
					var.writeFolder = value;
					break;		
				case "correcttotalreadcount":
					var.correctTotalReadCount = Double.parseDouble(value);
					break;
				case "spearman":
					var.spearman = Double.parseDouble(value);
					break;
				case "rlog":
					var.rLog = Double.parseDouble(value);
					break;
				case "addbeforelog":
					var.addLogVal = Double.parseDouble(value);
					break;
				case "lowesttoaverage":
					var.setLowestToAverage = Boolean.parseBoolean(value);
					break;
				case "adjustsampleaverages":
					var.adjustSampleAverages = Boolean.parseBoolean(value);
					break;
				default:
					System.out.println("Incorrect argument supplied:\n"+ args[a] +"\nexiting");
					System.exit(1);
			}
		}
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