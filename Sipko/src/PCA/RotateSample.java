package PCA;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;

import pca.MatrixStruct;
import pca.PCA;

public class RotateSample {

	public static void main(String[] args) throws IOException 
	{
//		String sampleFile = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/" + "RandomSamples.txt";
//		String vectorFolder = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/TESTexpression/";
		
//		String sampleFile = "E:/Groningen/Data/PublicSamples/Test7/" + "4DownSyndrome3Normal3Cancer_countsRND.txt";
//		String vectorFolder = "E:/Groningen/Data/PublicSamples/Test7/est_counts_nocancernocelllineSTdevCorrected/";
//		
//		boolean tpm = false;
//		double correctTotalReadCount = 0;
//		
//		checkArgs(args);
//		if(!System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko") || args.length !=0)
//		{
//			sampleFile = args[0];
//			vectorFolder = args[1]+"/";
//		}
//		String writeFolder = sampleFile.replace(".txt", "_Adj/");
//		
//		rotate(sampleFile, vectorFolder, writeFolder, tpm, correctTotalReadCount);	
	}
//	public static MatrixStruct[] rotate(String sampleFile, String vectorFolder, String writeFolder, boolean tpm, double correctTotalReadCount) throws IOException
//	{
//		return rotate(sampleFile, vectorFolder, writeFolder, true, true, tpm, correctTotalReadCount, -1, true);
//	}
	
	public static MatrixStruct[] rotate(String sampleFile, String vectorFolder, String writeFolder, boolean STdevCorrect,
			boolean log2, boolean tpm, double correctTotalReadCount, double spearman, boolean adjustSampleAverages,
			boolean setLowestToAverage, double addBeforeLog, double rLog, String singleSampleFN, String chromLocationsFile) throws IOException
	{
		/**6. Calculate PCscores for single sample**/
		pca.PCA.log(" 1. Loading sample matrix");
		MatrixStruct singleSample = new MatrixStruct(sampleFile);//expressionMatrix.getRow(0);
		singleSample.putGenesOnRows();
		
		String columnAveragesFN = vectorFolder+"SAMPLE_Norm_columnAverages.txt";
		MatrixStruct columnAverages = new MatrixStruct(columnAveragesFN);
		singleSample = center(singleSample, vectorFolder, writeFolder, STdevCorrect, log2, tpm, 
				correctTotalReadCount,spearman, columnAverages, adjustSampleAverages, setLowestToAverage,
				addBeforeLog, rLog, singleSampleFN);
		String saveNameSingleSampleScore = writeFolder + "SAMPLE.PC.txt";
		pca.PCA.log(" 11. Loading gene eigen vector file: ");
		MatrixStruct geneEigenVectors = new MatrixStruct(vectorFolder+"GENE.eigenvectors.txt", -1, 1000);//maximum 1000 PCs
		if(geneEigenVectors.getColHeaders()[0].contains("ENSG") || geneEigenVectors.getColHeaders()[0].contains("ENST"))
		{
			pca.PCA.log("Transposing");
			geneEigenVectors.transpose();
		}
		
		columnAverages.keepRows(geneEigenVectors);//also changes the rows of geneEigenvectors to have the same positions as columnAverages
		System.out.println("sample " + singleSample.rows() + " eigenvectors = " + geneEigenVectors.rows());
		singleSample.putGenesOnRows();
		singleSample.keepRows1Matrix(geneEigenVectors);//necessary for some normalizations
		System.out.println("sample1 " + singleSample.rows() + " eigenvectors = " + geneEigenVectors.rows());

		pca.PCA.log("Transposing");
		geneEigenVectors.transpose();
		geneEigenVectors.write(vectorFolder+"GENE.eigenvectors2.txt");

		pca.PCA.log(" 12. Rotating file to PC space: ");
		MatrixStruct[] scoreResults = PCA.scores(geneEigenVectors,singleSample, saveNameSingleSampleScore,false, false);
		
		pca.PCA.log("Files Written to: " + writeFolder);
		return scoreResults;
	}
//	public static MatrixStruct center(MatrixStruct singleSample, String vectorFolder, String writeFolder, MatrixStruct columnAverages) throws IOException
//	{
//		return center(singleSample, vectorFolder, writeFolder, true, true,false,0, -1, columnAverages, true);
//	}
	public static MatrixStruct center(MatrixStruct singleSample, String vectorFolder, String writeFolder, 
			boolean STdevCorrect, boolean log2,boolean tpm, double correctTotalReadCount, double spearman, 
			MatrixStruct columnAverages, boolean adjustSampleAverages, boolean setLowestToAverage,
			double addBeforeLog, double rLog, String singleSampleFN) throws IOException 
	{
		makeFolder(writeFolder);
		pca.PCA.log(" 3. Removing rows that do not exist in the averages vector from public data");
		
		if(!tpm && correctTotalReadCount < 1)
		{
			MatrixStruct quantVector = new MatrixStruct(vectorFolder+"SAMPLE_QuantileVector.txt");
			pca.PCA.log(" 4. Quantile normalization adjustion");
			singleSample.expressionToRank(quantVector);
			
			String quantileAdjustedFN = writeFolder+ "Quantile_adjusted.txt";	
			pca.PCA.log(" 5. Writing quantile normalization adjusted file to:" + quantileAdjustedFN);
			singleSample.write(quantileAdjustedFN);
		}
		
		if(correctTotalReadCount >0)
		{
			pca.PCA.log(" 6. Correcting for total read count");
			//String correctedNotLogged =  writeFolder+ "SAMPLE_TotalReadCountNormalized.txt";
			singleSample.correctForTotalReadCount(correctTotalReadCount,0.5);
			singleSample.write(writeFolder + "correctTotalReadCount_"+correctTotalReadCount+".txt");
		}
		if(rLog > 0)
		{
			pca.PCA.log(" 6. rLog transformation");
			MatrixStruct geoMeans = new MatrixStruct(vectorFolder+"geoMean.txt");
			String swapFN = writeFolder + "swapFile.txt";
			singleSample.write(swapFN);
			singleSample.rLog(rLog, writeFolder, swapFN, geoMeans);
			singleSample.write(writeFolder + "rLogTransformed_"+rLog+".txt");
		}
		
		if(log2)
		{
			if(correctTotalReadCount <= 0 && rLog <= 0) // need to add 1 before log to avoid log(0)
				addBeforeLog = 0.5;
			pca.PCA.log(" 6. Log transforming");
			singleSample.log2Transform(addBeforeLog);//Doing this after the quantile normalization now
			singleSample.write(writeFolder + "normalized_log2.txt");
		}
		
		columnAverages.keepRows(singleSample);
		singleSample.write(writeFolder+"keepRows.txt");
		
		if(spearman >= 0)
		{
			pca.PCA.log("Changing expression data into rank data");
			expressionToPublicRank(singleSample,spearman,vectorFolder+ "beforeRanks.txt",vectorFolder+ "ranks.txt");
			singleSample.write(writeFolder+"RankValues.txt");
		}
		singleSample.transpose();
		
		if(STdevCorrect)
		{
			pca.PCA.log(" 7. Adjusting for STdevs");
			String stdevFile = vectorFolder+ "gene_STDevs.txt";
			MatrixStruct stDevs = new MatrixStruct(stdevFile);
			singleSample.divideBy(stDevs, false);
			singleSample.write(writeFolder+"Centered.DivideBySTdev.txt");
		}
		if(setLowestToAverage)//this will cause the lowest values not to contribute to correlation or covariance
			LowestToAverage.lowestToAverage(singleSample);
		
		//System.out.println("singleSample=" + singleSample.matrix.get(1, 1));
		pca.PCA.log(" 8. Adjusting for column averages (centering to target PC space)");
		singleSample.adjustForAverageAllCols(columnAverages);
		//System.out.println("singleSample=" + singleSample.matrix.get(1, 1));
		String centeredColsFN = writeFolder+ "centeredColsOnly.txt";
		singleSample.write(centeredColsFN);
		
		MatrixStruct averages = singleSample.getAveragesPerRow();
		String rowAveragesFileName = writeFolder+"rowAverages.txt";
		averages.write(rowAveragesFileName);
		if(adjustSampleAverages)
		{
			pca.PCA.log(" 9. Adjusting for row averages (centering to target PC space)");
			singleSample.adjustForAverageAllrows(averages);
		}
		String centeredFN = writeFolder+ "centered.txt";
		
		pca.PCA.log(" 10. Writing PC centered file to: " + centeredFN);
		singleSample.write(centeredFN);
		return singleSample;
	}
	private static void expressionToPublicRank(MatrixStruct singleSample, double spearman, String beforeFN,
			String afterFN) {
		//
		MatrixStruct expression = new MatrixStruct(beforeFN);
		MatrixStruct rank = new MatrixStruct(afterFN);
		expression.transpose();//Put gene names on the rows
		expression.keepRows(singleSample);
		rank.transpose();
		rank.keepRows(singleSample);
		for(int c = 0; c < singleSample.cols(); c++)
		{
			for(int r = 0; r < singleSample.rows(); r++)
			{
				int index = getSmallestIndex(expression, singleSample,r,c);
				double rankValue = rank.matrix.get(r, index);
				singleSample.matrix.set(r, c, rankValue);
			}
		}
		
	}
	private static int getSmallestIndex(MatrixStruct expression, MatrixStruct singleSample, int r, int c) {
		double smallest = Double.POSITIVE_INFINITY;
		int smallestIndex = 0;
		for(int cRank = 0; cRank < expression.cols(); cRank++)
		{
			double publicVal = expression.matrix.get(r, cRank);
			double sampleVal = singleSample.matrix.get(r, c);
			double difference = Math.abs(publicVal-sampleVal);
			if(difference < smallest)
			{
				smallest = difference;
				smallestIndex = cRank;
			}
		}
		return smallestIndex;
	}
	public static void checkArgs(String[] args)
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
			return;
		if(args.length!= 2)
		{
			System.out.println("This script rotates the data into a supplied PC space (Defined by argument 2):\n"
					+ "It uses the following 2 arguments:\n"
					+ "1. expression File of the samples (samples on row X genes on columns) \n"
					+ "2. Folder wehre the vector and corresponding files, required for the rotation ar elocated\n"
					+ "   This folder can be created using CreateGeneEigenvectorFile.java script\n");
			System.exit(1);
		}
	}
	static void makeFolder(String writeFolder) 
	{
		File folder = new File(writeFolder);
		if(!folder.exists())
		{
			folder.mkdir();
		}
		
	}
}
