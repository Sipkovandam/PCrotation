//package PCA;
//
//import java.io.File;
//import java.io.IOException;
//import java.nio.file.Paths;
//import java.util.Hashtable;
//
//import JuhaPCA.PCA;
//import MatrixScripts.MyMatrix;
//import MatrixScripts.Zscore;
//import eqtlmappingpipeline.normalization.Normalizer;
//import umcg.genetica.math.matrix.DoubleMatrixDataset;
//
//public class RotateSample
//{
//
//	//"Rotates" a sample into the same space as another sample, by normalizing and centering it in the same way.
//
//	public static void main(String[] args) throws IOException
//	{
//		//		String sampleFile = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/" + "RandomSamples.txt";
//		//		String vectorFolder = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/TESTexpression/";
//
//		//		String sampleFile = "E:/Groningen/Data/PublicSamples/Test7/" + "4DownSyndrome3Normal3Cancer_countsRND.txt";
//		//		String vectorFolder = "E:/Groningen/Data/PublicSamples/Test7/est_counts_nocancernocelllineSTdevCorrected/";
//		//		
//		//		boolean tpm = false;
//		//		double correctTotalReadCount = 0;
//		//		
//		//		checkArgs(args);
//		//		if(!System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko") || args.length !=0)
//		//		{
//		//			sampleFile = args[0];
//		//			vectorFolder = args[1]+"/";
//		//		}
//		//		String writeFolder = sampleFile.replace(".txt", "_Adj/");
//		//		
//		//		rotate(sampleFile, vectorFolder, writeFolder, tpm, correctTotalReadCount);	
//	}
//	public static MyMatrix[] rotate(PcaPipeline var) throws IOException
//	{
//		return rotate(var, 1002);
//	}
//	public static MyMatrix[] rotate(PcaPipeline var, int maxPCs) throws IOException
//	{
//		/** 6. Calculate PCscores for single sample **/
//		JuhaPCA.PCA.log(" 1. Loading sample matrix");
//		MyMatrix singleSample = new MyMatrix(var.sampleFile);//expressionMatrix.getRow(0);
//		singleSample.putGenesOnRows();
//
//		//keep only the genes/rows that were used in the public samples as well
//		String geneAveragesFN = var.writeFolder + "SAMPLE_Norm_GeneAverages.txt";
//		MyMatrix geneAverages = new MyMatrix(geneAveragesFN);
//		geneAverages.keepRows(singleSample);
//
//		//rotate the sample to the same sample space
//		singleSample = center(	singleSample,
//								geneAverages,
//								var);
//
//		//calculate the PC scores
//		//Get the eigenvector file based on the public data and put it in the right orientation.
//		String saveNameSingleSampleScore = var.writeFolderCorrected + "SAMPLE.PC.txt";
//		JuhaPCA.PCA.log(" 11. Loading gene eigen vector file: ");
//		
//		File eigenVectorFn =null;
//		eigenVectorFn = new File(var.writeFolder + "eigenvectors.txt.gz");
//		if (!eigenVectorFn.exists())
//			eigenVectorFn = new File(var.writeFolder + "eigenvectors.txt");
//		
//
//		//read in matrix (genes are on the rows, PCs are on the columns in the eigenFile)
//		MyMatrix eigenVectors = new MyMatrix(	eigenVectorFn.getAbsolutePath(),
//													-1,
//													maxPCs+1);
//
//		if (eigenVectors.getRowHeaders()[0].contentEquals("PC1"))
//			eigenVectors.transpose();
//
//		if(var.isPcaOverGenes)
//			geneAverages.keepRows(eigenVectors);//also changes the rows of geneEigenvectors to have the same positions as geneAverages
//		
//		if(var.isPcaOverGenes)
//			singleSample.putGenesOnRows();
//		else if(!hasAnyOverlap(singleSample.getRowHash(), eigenVectors.rowNames))
//			singleSample.transpose();
//			
//		if(!hasAnyOverlap(singleSample.getRowHash(), eigenVectors.rowNames))
//		{
//			System.out.println("Error: No resemblant rownames between eigenvector and sample rownames");
//			System.exit(2);
//		}
//		
//		singleSample.keepRows1Matrix(eigenVectors);//necessary with some normalizations (quantile norm can introduce more 0 values into your matrix after the normalization creating a need for these to be removed afterward again)
//		
//		if (!eigenVectors.getRowHeaders()[0].contentEquals("PC1"))
//			eigenVectors.transpose();
//		System.out.println("geneEigenVectors.rows() = " + eigenVectors.rows() + " cols = " + eigenVectors.cols());
//		JuhaPCA.PCA.log(" 12. Calculate the PC scores: ");
//		MyMatrix[] scoreResults = PCA.scores(	eigenVectors,
//													singleSample,
//													saveNameSingleSampleScore,
//													false,
//													false);
//
//		JuhaPCA.PCA.log("Files Written to: " + var.writeFolderCorrected);
//		return scoreResults;
//	}
//
//	private static boolean hasAnyOverlap(	Hashtable<String, Integer> inputHash,
//											String[] names)
//	{
//		for(String name: names)
//		{
//			if(inputHash.containsKey(name))
//				return true;
//		}
//		return false;
//	}
//
//	//	public static Matrix center(Matrix singleSample, String vectorFolder, String writeFolder, Matrix columnAverages) throws IOException
//	//	{
//	//		return center(singleSample, vectorFolder, writeFolder, true, true,false,0, -1, columnAverages, true);
//	//	}
//	public static MyMatrix center(	MyMatrix singleSample,
//										MyMatrix geneAverages,
//										PcaPipeline var) throws IOException
//	{
//		makeFolder(var.writeFolderCorrected);
//		JuhaPCA.PCA.log(" 3. Removing rows that do not exist in the averages vector from public data");
//		geneAverages.keepRows(singleSample);
//		if (var.quantileNorm && var.correctTotalReadCount < 1)//prevents quantile norm from happening whilst also correctTotalReadCountis happening
//		{
//			MyMatrix quantVector = new MyMatrix(var.writeFolder + "SAMPLE_QuantileVector.txt");
//			JuhaPCA.PCA.log(" 4. Quantile normalization adjustion");
//			//use quantile normalize distribtion from public data <quantVector>
//			singleSample.expressionToRank(quantVector);
//
//			String quantileAdjustedFN = var.writeFolderCorrected + "Quantile_adjusted.txt";
//			JuhaPCA.PCA.log(" 5. Writing quantile normalization adjusted file to:" + quantileAdjustedFN);
//			singleSample.write(quantileAdjustedFN);
//		}
//
//		//correct for the total number of reads in a sample
//		if (var.correctTotalReadCount > 0)
//		{
//			JuhaPCA.PCA.log(" 6. Correcting for total read count");
//			//String correctedNotLogged =  writeFolder+ "SAMPLE_TotalReadCountNormalized.txt";
//			singleSample.correctForTotalReadCount(	var.correctTotalReadCount,
//													0.5);
//			singleSample.write(var.writeFolderCorrected + "correctTotalReadCount_" + var.correctTotalReadCount + ".txt");
//		}
//		if (var.deSeqNorm)
//		{
//			JuhaPCA.PCA.log(" 6. rLog transformation");
//			MyMatrix geoMeans = new MyMatrix(var.writeFolder + "geoMean.txt");
//			String swapFN = var.writeFolderCorrected + "swapFile.txt";
//			singleSample.write(swapFN);
//			DeSeqNorm.rLog(	singleSample,
//						var.writeFolderCorrected,
//						swapFN,
//						geoMeans,
//						null);
//			singleSample.write(var.writeFolderCorrected + "rLogTransformed_" + var.deSeqNorm + ".txt");
//		}
//
//		if (var.log2)
//		{
//			//			if(correctTotalReadCount <= 0 && rLog <= 0) // need to add 1 before log to avoid log(0)
//			//				addLogVal = 0.5;
//			JuhaPCA.PCA.log(" 6. Log transforming");
//			singleSample.log2Transform(var.addLogVal);//Doing this after the quantile normalization now
//			singleSample.write(var.writeFolderCorrected + "normalized_log2.txt");
//		}
//
//		singleSample.write(var.writeFolderCorrected + "keepRows.txt");
//
//		if (var.spearman >= 0)
//		{
//			JuhaPCA.PCA.log("Changing expression data into rank data");
//			expressionToPublicRank(	var,
//									singleSample,
//									var.spearman,
//									var.writeFolder + "beforeRanks.txt",
//									var.writeFolder + "ranks.txt");
//			singleSample.write(var.writeFolderCorrected + "RankValues.txt");
//		}
//
//		if (var.correctInputForSTdevs)
//		{
//			JuhaPCA.PCA.log(" 7. Adjusting for STdevs");
//			String stdevFile = var.writeFolder + "gene_STDevs.txt";
//			MyMatrix stDevs = new MyMatrix(stdevFile);
//			singleSample.divideBy(	stDevs,
//									false);
//			singleSample.write(var.writeFolderCorrected + "Centered.DivideBySTdev.txt");
//		}
//		if (var.setLowestToAverage)//this will cause the lowest values not to contribute to correlation or covariance
//			LowestToAverage.lowestToAverage(singleSample);
//
//		//calculate z-scores for the input matrix based on avgStdevFolder+"avgStdev_uncorrected.txt" file. If avgStdevFolder==null, calculates avgs and stdevs per gene based on current matrix.
//		if(var.zScores)
//			Zscore.zScores(	var.writeFolderCorrected,
//							"input_zScores",
//							singleSample,
//							var.avgStdevFolder,
//							"avgStdev_normLog2.txt",
//							false);
//		if (var.centerGenes)
//		{
//			JuhaPCA.PCA.log(" 9. Adjusting for gene averages (centering to target PC space)");
//			singleSample.adjustForAverageAllGenes(geneAverages);
//		}
//		
//		MyMatrix sampleAvgs = singleSample.getAveragesPerCol();
//		String sampleAveragesFileName = var.writeFolderCorrected + "sampleAverages.txt";
//		sampleAvgs.write(sampleAveragesFileName);
//		if (var.centerSamples)
//		{
//			JuhaPCA.PCA.log(" 8. Adjusting for row averages (centering to target PC space)");
//			singleSample.adjustForAverageAllSamples(sampleAvgs);
//		}
//
//		String centeredFN = var.writeFolderCorrected + "centered.txt";
//		JuhaPCA.PCA.log(" 10. Writing PC centered file to: " + centeredFN);
//		singleSample.write(centeredFN);
//		return singleSample;
//	}
//
//	private static void expressionToPublicRank(	PcaPipeline var,
//												MyMatrix singleSample,
//												double spearman,
//												String beforeFN,
//												String afterFN)
//	{
//		//
//		MyMatrix expression = new MyMatrix(beforeFN);
//		MyMatrix rank = new MyMatrix(afterFN);
//		expression.putGenesOnRows();//Put gene names on the rows
//		expression.keepRows(singleSample);
//		rank.transpose();
//		rank.keepRows(singleSample);
//		for (int c = 0; c < singleSample.cols(); c++)
//		{
//			for (int r = 0; r < singleSample.rows(); r++)
//			{
//				int index = getSmallestIndex(	expression,
//												singleSample,
//												r,
//												c);
//				double rankValue = rank.matrix.get(	r,
//													index);
//				singleSample.matrix.set(r,
//										c,
//										rankValue);
//			}
//		}
//
//	}
//
//	private static int getSmallestIndex(MyMatrix expression,
//										MyMatrix singleSample,
//										int r,
//										int c)
//	{
//		double smallest = Double.POSITIVE_INFINITY;
//		int smallestIndex = 0;
//		for (int cRank = 0; cRank < expression.cols(); cRank++)
//		{
//			double publicVal = expression.matrix.get(	r,
//														cRank);
//			double sampleVal = singleSample.matrix.get(	r,
//														c);
//			double difference = Math.abs(publicVal - sampleVal);
//			if (difference < smallest)
//			{
//				smallest = difference;
//				smallestIndex = cRank;
//			}
//		}
//		return smallestIndex;
//	}
//
//	public static void checkArgs(String[] args)
//	{
//		if (System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
//			return;
//		if (args.length != 2)
//		{
//			System.out.println("This script rotates the data into a supplied PC space (Defined by argument 2):\n" + "It uses the following 2 arguments:\n" + "1. expression File of the samples (samples on row X genes on columns) \n" + "2. Folder wehre the vector and corresponding files, required for the rotation ar elocated\n" + "   This folder can be created using CreateGeneEigenvectorFile.java script\n");
//			System.exit(1);
//		}
//	}
//
//	static void makeFolder(String writeFolder)
//	{
//		File folder = new File(writeFolder);
//		if (!folder.exists())
//		{
//			folder.mkdir();
//		}
//
//	}
//}
