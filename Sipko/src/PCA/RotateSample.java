package PCA;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;

import pca.MatrixStruct;
import pca.PCA;

public class RotateSample {

	public static void main(String[] args) throws IOException 
	{
		String sampleFile = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/" + "TESTexpression.txt";
		String vectorFolder = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/TESTexpression/";
		
		checkArgs(args);
		if(!System.getProperty("user.dir").contains("E:\\Groningen\\Workspace"))
		{
			sampleFile = args[0];
			vectorFolder = args[1]+"/";
		}
		String writeFolder = sampleFile.replace(".txt", "_Adj/");
		
		rotate(sampleFile, vectorFolder, writeFolder);	
	}
	public static MatrixStruct[] rotate(String sampleFile, String vectorFolder, String writeFolder) throws IOException
	{
		/**6. Calculate PCscores for single sample**/
		makeFolder(writeFolder);
		pca.PCA.log(" 1. Loading sample matrix");
		Matrix singleSample = new Matrix(sampleFile);//expressionMatrix.getRow(0);
		
		pca.PCA.log(" 2. Transposing");
		singleSample.transpose();
		Matrix quantVector = new Matrix(vectorFolder+"SAMPLE_QuantileVector.txt");
		pca.PCA.log(" 3. Removing rows that do not exist in the quantile normalization vector");
		singleSample.keepRows(quantVector);
		pca.PCA.log(" 4. Quantile normalization adjustion");
		singleSample.quantileNormAdjust(quantVector);
		String quantileAdjustedFN = writeFolder+ "Quantile_adjusted.txt";
		pca.PCA.log(" 5. Writing quantile normalization adjusted file to:" + quantileAdjustedFN);
		pca.PCA.log("Log transforming");
		singleSample.log2Transform();//Doing this after the quantile normalization now
		singleSample.write(quantileAdjustedFN);
		singleSample.transpose();
		pca.PCA.log(" 6. Adjusting for column averages (centering to target PC space)");
		singleSample.adjustForAverageAllCols(new Matrix(vectorFolder+"SAMPLE_QuantNorm_columnAverages.txt"));
		pca.PCA.log(" 7. Adjusting for row averages (centering to target PC space)");
		singleSample.adjustForAverageAllrows(new Matrix(vectorFolder+"SAMPLE_QuantNorm_rowAverages.txt"));
		String centeredFN = writeFolder+ "Quantile_adjusted.centered.txt";
		pca.PCA.log(" 8. Writing PC centered file to: " + centeredFN);
		singleSample.write(centeredFN);
		pca.PCA.log(" 9. Rotating file to PC space: ");
		MatrixStruct singleSampleStruct = new MatrixStruct(singleSample.rowNames, singleSample.colNames, singleSample.values);
		String averagesFN = vectorFolder + "GENE_PC.averages.txt";//used as input
		String saveNameSingleSampleScore = writeFolder + "SAMPLE.PC.txt";
		MatrixStruct geneEigenVectors = new MatrixStruct(vectorFolder+"GENE.eigenvectors.txt");
		MatrixStruct[] scoreResults = PCA.scores(geneEigenVectors,singleSampleStruct, Paths.get(averagesFN), saveNameSingleSampleScore);
		pca.PCA.log("Files Written to: " + writeFolder);
		return scoreResults;
	}
	public static void checkArgs(String[] args)
	{
		if(System.getProperty("user.dir").contains("E:\\Groningen\\Workspace"))
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
	private static void makeFolder(String writeFolder) 
	{
		File folder = new File(writeFolder);
		if(!folder.exists())
		{
			folder.mkdir();
		}
		
	}
}
