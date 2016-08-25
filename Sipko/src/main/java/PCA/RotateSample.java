package PCA;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;

import JuhaPCA.PCA;
import eqtlmappingpipeline.normalization.Normalizer;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

public class RotateSample {

	//"Rotates" a sample into the same space as another sample, by normalizing and centering it in the same way.
	
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
	
	public static MatrixStruct[] rotate(String sampleFile, 
										String vectorFolder, 
										String writeFolder, 
										boolean STdevCorrect,
										boolean log2,
										boolean skipQuantileNorm, 
										double correctTotalReadCount, 
										double spearman, 
										boolean adjustSampleAverages,
										boolean setLowestToAverage, 
										double addLogVal, 
										boolean rLog, 
										String singleSampleFN, 
										String chromLocationsFile, 
										String correctGC, Var var) throws IOException
	{
		/**6. Calculate PCscores for single sample**/
		JuhaPCA.PCA.log(" 1. Loading sample matrix");
		MatrixStruct singleSample = new MatrixStruct(sampleFile);//expressionMatrix.getRow(0);
		singleSample.putGenesOnRows();
		
		//keep only the genes/rows that were used in the public samples as well
		String geneAveragesFN = vectorFolder+"SAMPLE_Norm_GeneAverages.txt";
		MatrixStruct columnAverages = new MatrixStruct(geneAveragesFN);
		columnAverages.keepRows(singleSample);
		
		//rotate the sample to the same sample space
		singleSample = center(	singleSample, 
								vectorFolder, 
								writeFolder, 
								STdevCorrect, 
								log2, 
								skipQuantileNorm,
								correctTotalReadCount,
								spearman, 
								columnAverages, 
								adjustSampleAverages, 
								setLowestToAverage,
								addLogVal, 
								rLog, 
								singleSampleFN, 
								correctGC);
		
		//calculate the PC scores
		//Get the eigenvector file based on the public data and put it in the right orientation.
		String saveNameSingleSampleScore = writeFolder + "SAMPLE.PC.txt";
		JuhaPCA.PCA.log(" 11. Loading gene eigen vector file: ");
		File eigenFile = new File(vectorFolder+"GENE.eigenvectors.txt.gz");
		if(!eigenFile.exists())
			eigenFile = new File(vectorFolder+"GENE.eigenvectors.txt");
	
		MatrixStruct geneEigenVectors = null;
		String eigen2Name = vectorFolder+"GENE.eigenvectors2.txt.gz";
		if(!var.reUseEigen2 || !new File(eigen2Name).exists())
		{
			//read in matrix (genes are on the rows, PCs are on the columns in the eigenFile)
			geneEigenVectors = new MatrixStruct(eigenFile.getAbsolutePath(), -1, 5001);//maximum 5001 PCs
			
			//geneEigenVectors.putGenesOnRows();

		}
		else
		{
			geneEigenVectors = new MatrixStruct(eigen2Name);//maximum 5001 PCs
		}
		System.out.println("sample0 " + singleSample.rows() + " eigenvectors = " + geneEigenVectors.rows());
		columnAverages.keepRows(geneEigenVectors);//also changes the rows of geneEigenvectors to have the same positions as columnAverages
		System.out.println("sample " + singleSample.rows() + " eigenvectors = " + geneEigenVectors.rows());
		singleSample.putGenesOnRows();
		singleSample.keepRows1Matrix(geneEigenVectors);//necessary with some normalizations
		System.out.println("sample1 " + singleSample.rows() + " eigenvectors = " + geneEigenVectors.rows());
		
		//if the eigen2 file does not exist, write it so we can do this faster next time we correct our data
		if(!new File(eigen2Name).exists())
		{
			JuhaPCA.PCA.log("Transposing");
			geneEigenVectors.transpose();//Puts the genes on the columns (and the PCs on the rows)
			geneEigenVectors.write(vectorFolder+"GENE.eigenvectors2.txt");//had genes on columns and is only part of geneEigenvectors.txt (5000 PCs)
		}
		
		JuhaPCA.PCA.log(" 12. Calculate the PC scores: ");
		singleSample.write(vectorFolder+"TESTSAMPLEOUT.txt");
		geneEigenVectors.write(vectorFolder+"TESTVECTOROUT.txt");
		MatrixStruct[] scoreResults = PCA.scores(geneEigenVectors,singleSample, saveNameSingleSampleScore,false, false);
		
		JuhaPCA.PCA.log("Files Written to: " + writeFolder);
		return scoreResults;
	}
//	public static MatrixStruct center(MatrixStruct singleSample, String vectorFolder, String writeFolder, MatrixStruct columnAverages) throws IOException
//	{
//		return center(singleSample, vectorFolder, writeFolder, true, true,false,0, -1, columnAverages, true);
//	}
	public static MatrixStruct center(	MatrixStruct singleSample, 
										String vectorFolder, 
										String writeFolder, 
										boolean STdevCorrect, 
										boolean log2,
										boolean skipQuantileNorm, 
										double correctTotalReadCount, 
										double spearman, 
										MatrixStruct geneAverages, 
										boolean adjustSampleAverages,
										boolean setLowestToAverage,
										double addLogVal, 
										boolean rLog, 
										String singleSampleFN, 
										String correctGC) throws IOException 
	{
		makeFolder(writeFolder);
		JuhaPCA.PCA.log(" 3. Removing rows that do not exist in the averages vector from public data");
		geneAverages.keepRows(singleSample);
		if(!skipQuantileNorm && correctTotalReadCount < 1)//prevents quantile norm from happening whilst also correctTotalReadCountis happening
		{
			
			MatrixStruct quantVector = new MatrixStruct(vectorFolder+"SAMPLE_QuantileVector.txt");
			JuhaPCA.PCA.log(" 4. Quantile normalization adjustion");
			//use quantile normalize distribtion from public data <quantVector>
			singleSample.expressionToRank(quantVector);
			
			String quantileAdjustedFN = writeFolder+ "Quantile_adjusted.txt";	
			JuhaPCA.PCA.log(" 5. Writing quantile normalization adjusted file to:" + quantileAdjustedFN);
			singleSample.write(quantileAdjustedFN);
		}
		
		//correct for the total number of reads in a sample
		if(correctTotalReadCount >0)
		{
			JuhaPCA.PCA.log(" 6. Correcting for total read count");
			//String correctedNotLogged =  writeFolder+ "SAMPLE_TotalReadCountNormalized.txt";
			singleSample.correctForTotalReadCount(correctTotalReadCount,0.5);
			singleSample.write(writeFolder + "correctTotalReadCount_"+correctTotalReadCount+".txt");
		}
		if(rLog)
		{
			JuhaPCA.PCA.log(" 6. rLog transformation");
			MatrixStruct geoMeans = new MatrixStruct(vectorFolder+"geoMean.txt");
			String swapFN = writeFolder + "swapFile.txt";
			singleSample.write(swapFN);
			RLog.rLog(singleSample, writeFolder, swapFN, geoMeans,null);
			singleSample.write(writeFolder + "rLogTransformed_"+rLog+".txt");
		}
		
		if(log2)
		{
//			if(correctTotalReadCount <= 0 && rLog <= 0) // need to add 1 before log to avoid log(0)
//				addLogVal = 0.5;
			JuhaPCA.PCA.log(" 6. Log transforming");
			singleSample.log2Transform(addLogVal);//Doing this after the quantile normalization now
			singleSample.write(writeFolder + "normalized_log2.txt");
		}
		
		singleSample.write(writeFolder+"keepRows.txt");
		
		if(spearman >= 0)
		{
			JuhaPCA.PCA.log("Changing expression data into rank data");
			expressionToPublicRank(singleSample,spearman,vectorFolder+ "beforeRanks.txt",vectorFolder+ "ranks.txt");
			singleSample.write(writeFolder+"RankValues.txt");
		}
		
		if(STdevCorrect)
		{
			JuhaPCA.PCA.log(" 7. Adjusting for STdevs");
			String stdevFile = vectorFolder+ "gene_STDevs.txt";
			MatrixStruct stDevs = new MatrixStruct(stdevFile);
			singleSample.divideBy(stDevs, false);
			singleSample.write(writeFolder+"Centered.DivideBySTdev.txt");
		}
		if(setLowestToAverage)//this will cause the lowest values not to contribute to correlation or covariance
			LowestToAverage.lowestToAverage(singleSample);
		
		MatrixStruct averages = singleSample.getAveragesPerRow();
		String rowAveragesFileName = writeFolder+"rowAverages.txt";
		averages.write(rowAveragesFileName);
		
		//correct for GC content
		if(correctGC != null)
		{
			MatrixStruct gcPerGene = new MatrixStruct(correctGC);
			singleSample = GCcontent.calculateAndCorrect(singleSample,gcPerGene, writeFolder+"gCperSampleWriteFN.txt", writeFolder + "GCcorrected.txt.gz");
		}
		
		if(adjustSampleAverages)
		{
			JuhaPCA.PCA.log(" 8. Adjusting for row averages (centering to target PC space)");
			singleSample.adjustForAverageAllsamples(averages);
		}
		
		JuhaPCA.PCA.log(" 9. Adjusting for gene averages (centering to target PC space)");
		singleSample.adjustForAverageAllGenes(geneAverages);
		String centeredFN = writeFolder+ "centered.txt";
		
		JuhaPCA.PCA.log(" 10. Writing PC centered file to: " + centeredFN);
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
