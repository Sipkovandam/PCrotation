package PCA;

import java.io.File;
import umcg.genetica.math.stats.concurrent.ConcurrentCovariation;
import java.io.IOException;
import java.nio.file.Paths;

import no.uib.cipr.matrix.NotConvergedException;
import pca.MatrixStruct;
import pca.PCA;

public class CreateGeneEigenvectorFile 
{
	public static void main(String[] args) throws IOException, NotConvergedException
	{
		String expFile = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/" + "TESTexpression.txt";
		String chromLocationsFile = "E:/Groningen/Data/GenePositionInfo_23X_24Y_25MT_26rest.txt";
		//String expFile = "E:/Groningen/Data/PublicSamples/ComputerTest/Rtest/" + "rSample.txt";
		System.out.println(System.getProperty("user.dir"));
		checkArgs(args);
		
		boolean writeAll = true;
		
		if(!System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
		{
			expFile = args[0];
			if(args.length>1)
				chromLocationsFile = args[1];
			else
				chromLocationsFile = null;
			if(args.length>2 && args[2].contains("nowrite"))
				writeAll = false;
		}
		
		run(expFile,chromLocationsFile, writeAll);
	}
	private static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
			return;
		if(args.length < 1)
		{
			System.out.println("This script calculates the eigenvectors over the genes and uses the following input:\n"
					+ "1. Expression file (samples on rows, gene names on columns\n"
					+ "(Optional) 2. File with genes names in 1st column, chromosome in 2nd, position in 3rd");
			System.exit(1);
		}
	}
	public static void run(String expFile, String chromLocationsFile, boolean writeAll) throws IOException, NotConvergedException
	{
		String writeFolder = expFile.replace(".txt", "/");
		makeFolder(writeFolder);
		Matrix expressionMatrix = null;
		MatrixStruct covMat = null;
		String expNormLogCentFile= writeFolder+"MATRIX_Centered.txt";
		if(!new File(writeFolder+"MATRIX_Centered.CovariationMatrix.txt").exists())
		{
			if(chromLocationsFile != null)
			{
				pca.PCA.log(" 0. Sorting IDs based on chromosome locations");
				MatrixStruct expressionMatrixStruct = new MatrixStruct(expFile);
				expressionMatrixStruct.transpose();
				MatrixStruct chromLocs = new MatrixStruct(chromLocationsFile);
				expressionMatrixStruct = SortChromosome.sort(expressionMatrixStruct, chromLocs);
				expressionMatrix = new Matrix(expressionMatrixStruct);
				expressionMatrixStruct = null;
			}
			else
			{
				expressionMatrix = new Matrix(expFile);
				pca.PCA.log(" 1. Transposing");
				expressionMatrix.transpose();
			}
			expressionMatrix.write(writeFolder+"startMatrix.txt");
			pca.PCA.log(" 2. Removing genes without variance");
			expressionMatrix.removeNoVariance();
			String removedGenesFN = writeFolder+"noVarRemoved.txt";
			if(writeAll)expressionMatrix.write(removedGenesFN);
			pca.PCA.log(" 3. Calculating quantile normalization vector");
			Matrix qNormVector = expressionMatrix.quantileNormVector();
			if(writeAll)qNormVector.write(writeFolder+ "SAMPLE_QuantileVector.txt");
			pca.PCA.log(" 4. Quantile normalization");//i changed this to come after the Quantile normalization
			expressionMatrix.quantileNormAdjust(qNormVector);
			String quantFNnotLogged = writeFolder+ "SAMPLE_QuantileNormalized.txt";
			pca.PCA.log(" 5. Writing quantile normalized data in: " + quantFNnotLogged);
			if(writeAll)expressionMatrix.write(quantFNnotLogged);
			pca.PCA.log(" 6. Log2 transforming");
			expressionMatrix.log2Transform();
			String quantFN = writeFolder+ "SAMPLE_QuantileNormalizedLog2.txt";
			pca.PCA.log(" 7. Writing logged SAMPLE_QuantileNormalizedLog2 normalized data in: " + quantFN);
			if(writeAll)expressionMatrix.write(quantFN);
			pca.PCA.log(" 8. Transposing");
			expressionMatrix.transpose();
			//Centering is done before the PCA inside the PCA function, to assure the centering is done over the right dimension (PCA function also puts matrixes in right orientation)
			pca.PCA.log(" 9. Calculating column averages");
			Matrix colAverages = expressionMatrix.calcAvgCols();//NOTE: this is a vector that is based on rows already being corrected for averages.
			colAverages.write(writeFolder+ "SAMPLE_QuantNorm_columnAverages.txt");
			pca.PCA.log("10. Centering: Adjusting for column averages");
			expressionMatrix.adjustForAverageAllCols(colAverages);
			pca.PCA.log("10. Calculating row averages");
			Matrix rowAverages = expressionMatrix.calcAvgRows();
			rowAverages.write(writeFolder+ "SAMPLE_QuantNorm_rowAverages.txt");
			pca.PCA.log("11. Adjusting for row averages");
			expressionMatrix.adjustForAverageAllrows(rowAverages);
			//String expNormLogCentFile= writeFolder+"MATRIX_CenteredOverRowsAndColumns.txt";
			pca.PCA.log("11. Writing centered file in: " + expNormLogCentFile);
			expressionMatrix.write(expNormLogCentFile);
			
			
			pca.PCA.log("12. creating covariance matrix over the samples");
			//this function is slow ...
			covMat = pca.PCA.covariance(Paths.get(expNormLogCentFile), writeFolder+"SAMPLE_covariance.txt");
			covMat.write(writeFolder+"MATRIX_Centered.CovJuhaScript.txt");
			//Use Marc Jan's script instead...
		    //java  -Xmx30g -jar javaPCoA-1.0.4-SNAPSHOT-jar-with-dependencies.jar -m covariation -r -t 20 -i /Volumes/Promise_RAID/GeneNetwork/Sipko/est_counts_nocancernocellline/MATRIX_CenteredOverRowsAndColumns.txt
//			ConcurrentCovariation calculator = new ConcurrentCovariation(20);
//			double[][] covMatrix = calculator.pairwiseCovariation(expressionMatrix.values);
//			covMat = new MatrixStruct(expressionMatrix.rowNames ,expressionMatrix.rowNames ,covMatrix);
//			covMat.write(writeFolder+"SAMPLE_covariance.txt");
		}  	
		else
		{
			pca.PCA.log("13.1. Reading covariance matrix");
			covMat = new MatrixStruct(writeFolder+"MATRIX_Centered.CovariationMatrix.txt");
			pca.PCA.log("13.2. Reading adjusted expression matrix");
			expressionMatrix =  new Matrix(expNormLogCentFile);
			
		}
		
		pca.PCA.log("14. calculating eigenvalues over the samples");
		MatrixStruct[] evds = pca.PCA.evd(covMat, Paths.get(writeFolder+"SAMPLE"));
		MatrixStruct eigenVectors = evds[0];
		MatrixStruct PCeigenvalues = evds[1];
		
		/**6. Calculate PCscores of the genes of public/healthy expression data**/
		pca.PCA.log("15. calculating PCscores over the genes");
		MatrixStruct expressionStruct = new MatrixStruct(expressionMatrix.rowNames, expressionMatrix.colNames, expressionMatrix.values);
		String saveNamePCscoresGene = writeFolder + "GENE_PC.txt";
		MatrixStruct[] PCscoresGenesAndAverages = PCA.scores(eigenVectors,expressionStruct, null, saveNamePCscoresGene);
		MatrixStruct PCscoresGenes = PCscoresGenesAndAverages[0];
		//MatrixStruct averages = PCAscoresSamplesAndAverages[1];
		System.gc();
		
		/**7. Transform PCscores to eigenvectors of Genes**/
		pca.PCA.log("16. Transform PCscores to eigenvectors of Genes");
		String saveNameEigenVectorsOverGenes = writeFolder + "GENE.eigenvectors.txt";
		MatrixStruct geneEigenVectors = PCA.transform(PCscoresGenes, PCeigenvalues, saveNameEigenVectorsOverGenes);
		/**8. Calculating PCscores for all samples**/
		pca.PCA.log("17. Calculating PCscores for all samples");
		String averagesFN = writeFolder + "GENE_PC.averages.txt";//used as input
		MatrixStruct[] scoreResults = PCA.scores(geneEigenVectors,expressionStruct, Paths.get(averagesFN), writeFolder+"SAMPLE_PC.txt");
		MatrixStruct PCsampleScores = scoreResults[0];
		pca.PCA.log("18. Calculating Zscores for all PCscores for all samples");
		MatrixStruct zScoreStats = Zscore.changeToZscores(PCsampleScores);
		PCsampleScores.write(writeFolder+ "pcZscoresSamples.txt");
		zScoreStats.write(writeFolder+ "pcZscores_Stats.txt");
		
		pca.PCA.log("Files written to: " + writeFolder);
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
