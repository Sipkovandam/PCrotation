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
		
		String sampleFile = "E:/Groningen/Data/PublicSamples/Test7/" + "4DownSyndrome3Normal3Cancer_countsRND.txt";
		String vectorFolder = "E:/Groningen/Data/PublicSamples/Test7/est_counts_nocancernocelllineSTdevCorrected/";
		
		boolean tpm = false;
		boolean correctTotalReadCount = false;
		
		checkArgs(args);
		if(!System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko") || args.length !=0)
		{
			sampleFile = args[0];
			vectorFolder = args[1]+"/";
		}
		String writeFolder = sampleFile.replace(".txt", "_Adj/");
		
		rotate(sampleFile, vectorFolder, writeFolder, tpm, correctTotalReadCount);	
	}
	public static MatrixStruct[] rotate(String sampleFile, String vectorFolder, String writeFolder, boolean tpm, boolean correctTotalReadCount) throws IOException
	{
		return rotate(sampleFile, vectorFolder, writeFolder, true, true, tpm, correctTotalReadCount);
	}
	
	public static MatrixStruct[] rotate(String sampleFile, String vectorFolder, String writeFolder, boolean STdevCorrect, boolean log2, boolean tpm, boolean correctTotalReadCount) throws IOException
	{
		/**6. Calculate PCscores for single sample**/
		pca.PCA.log(" 1. Loading sample matrix");
		Matrix singleSample = new Matrix(sampleFile);//expressionMatrix.getRow(0);
		System.out.println("rows = " + singleSample.rowNames.length + " columns = " + singleSample.colNames.length);
		singleSample = center(singleSample, vectorFolder, writeFolder, STdevCorrect, log2, tpm, correctTotalReadCount);

		MatrixStruct singleSampleStruct = new MatrixStruct(singleSample.rowNames, singleSample.colNames, singleSample.values);
		String saveNameSingleSampleScore = writeFolder + "SAMPLE.PC.txt";
		pca.PCA.log(" 11. Loading gene eigen vector file: ");
		MatrixStruct geneEigenVectors = new MatrixStruct(vectorFolder+"GENE.eigenvectors.txt");
		pca.PCA.log(" 12. Rotating file to PC space: ");
		MatrixStruct[] scoreResults = PCA.scores(geneEigenVectors,singleSampleStruct, saveNameSingleSampleScore,false, false);
		
		pca.PCA.log("Files Written to: " + writeFolder);
		return scoreResults;
	}
	public static Matrix center(Matrix singleSample, String vectorFolder, String writeFolder) throws IOException
	{
		return center(singleSample, vectorFolder, writeFolder, true, true,false,false);
	}
	public static Matrix center(Matrix singleSample, String vectorFolder, String writeFolder, boolean STdevCorrect, boolean log2,boolean tpm, boolean correctTotalReadCount) throws IOException 
	{
		makeFolder(writeFolder);
		
		if(singleSample.colNames[0].contains("ENSG") || singleSample.colNames[0].contains("ENST"))
		{
			pca.PCA.log(" 2. Transposing");
			singleSample.transpose();
		}
		Matrix keepGenes = new Matrix(vectorFolder+"SAMPLE_Norm_columnAverages.txt");
		pca.PCA.log(" 3. Removing rows that do not exist in the quantile normalization vector");
		singleSample.keepRows(keepGenes);
		
		if(!tpm)
		{
			Matrix quantVector = new Matrix(vectorFolder+"SAMPLE_QuantileVector.txt");
			pca.PCA.log(" 4. Quantile normalization adjustion");
			singleSample.quantileNormAdjust(quantVector);
			
			String quantileAdjustedFN = writeFolder+ "Quantile_adjusted.txt";	
			pca.PCA.log(" 5. Writing quantile normalization adjusted file to:" + quantileAdjustedFN);
			singleSample.write(quantileAdjustedFN);
		}
		
		if(correctTotalReadCount)
		{
			pca.PCA.log(" 6. Correcting for total read count");
			String correctedNotLogged =  writeFolder+ "SAMPLE_TotalReadCountNormalized.txt";
			singleSample.correctForTotalReadCount();
		}
		
		if(log2)
		{
			pca.PCA.log(" 6. Log transforming");
			singleSample.log2Transform();//Doing this after the quantile normalization now
			singleSample.write(writeFolder + "normalized_log2.txt");
		}
		
		singleSample.transpose();

		if(STdevCorrect)
		{
			MatrixStruct singleSampleStruct = new MatrixStruct(singleSample.rowNames, singleSample.colNames, singleSample.values);
			pca.PCA.log(" 7. Adjusting for STdevs");
			String stdevFile = vectorFolder+ "gene_STDevs.txt";
			MatrixStruct stDevs = new MatrixStruct(stdevFile);
			singleSampleStruct.divideBy(stDevs, false);
			singleSampleStruct.write(writeFolder+"Centered.DivideBySTdev.txt");
			singleSample = new Matrix(singleSampleStruct);
		}
		
		pca.PCA.log(" 8. Adjusting for column averages (centering to target PC space)");
		singleSample.adjustForAverageAllCols(new Matrix(vectorFolder+"SAMPLE_Norm_columnAverages.txt"));
		String centeredColsFN = writeFolder+ "Quantile_adjusted.centeredColsOnly.txt";
		singleSample.write(centeredColsFN);
		
		pca.PCA.log(" 9. Adjusting for row averages (centering to target PC space)");
		Matrix averages = singleSample.calcAvgRows();
		String rowAveragesFileName = writeFolder+"rowAverages.txt";
		averages.write(rowAveragesFileName);
		singleSample.adjustForAverageAllrows(averages);
		String centeredFN = writeFolder+ "Quantile_adjusted.centered.txt";
		
		pca.PCA.log(" 10. Writing PC centered file to: " + centeredFN);
		singleSample.write(centeredFN);
		return singleSample;
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
