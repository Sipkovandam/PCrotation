package PCA;

import java.io.File;
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
		//String expFile = "E:/Groningen/Data/PublicSamples/ComputerTest/Rtest/" + "rSample.txt";
		System.out.println(System.getProperty("user.dir"));
		checkArgs(args);
		
		boolean writeAll = true;
		
		if(!System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
		{
			expFile = args[0];
			if(args.length ==2)
				writeAll = false;
		}
		
		run(expFile, writeAll);
	}
	private static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
			return;
		if(args.length!= 1)
		{
			System.out.println("This script calculates the eigenvectors over the genes and uses the following input:\n"
					+ "1. Expression file (samples on rows, gene names on columns");
			System.exit(1);
		}
	}
	public static void run(String expFile, boolean writeAll) throws IOException, NotConvergedException
	{
		String writeFolder = expFile.replace(".txt", "/");
		makeFolder(writeFolder);
		Matrix expressionMatrix = null;
		MatrixStruct covMat = null;
		String expNormLogCentFile= writeFolder+"MATRIX_CenteredOverRowsAndColumns.txt";
		if(!new File(writeFolder+"MATRIX_CenteredOverRowsAndColumns.CovariationMatrix.txt").exists())
		{
			expressionMatrix = new Matrix(expFile);
			pca.PCA.log(" 1. Transposing");
			expressionMatrix.transpose();
			pca.PCA.log(" 2. Removing genes without variance");
			expressionMatrix.removeNoVariance();
			String removedGenesFN = writeFolder+"noVarRemoved.txt";
			if(writeAll)expressionMatrix.write(removedGenesFN);
			pca.PCA.log(" 3. Calculating quantile normalization vector");
			Matrix qNormVector = expressionMatrix.quantileNormVector();
			if(writeAll)qNormVector.write(writeFolder+ "SAMPLE_QuantileVector.txt");
			pca.PCA.log(" 5. Quantile normalization");//i changed this to come after the Quantile normalization
			expressionMatrix.quantileNormAdjust(qNormVector);
			String quantFNnotLogged = writeFolder+ "SAMPLE_QuantileNormalized.txt";
			pca.PCA.log(" 6. Writing quantile normalized data in: " + quantFNnotLogged);
			if(writeAll)expressionMatrix.write(quantFNnotLogged);
			pca.PCA.log(" 4. Log2 transforming");
			expressionMatrix.log2Transform();
			String quantFN = writeFolder+ "SAMPLE_QuantileNormalizedLog2.txt";
			pca.PCA.log(" 6. Writing logged SAMPLE_QuantileNormalizedLog2 normalized data in: " + quantFN);
			if(writeAll)expressionMatrix.write(quantFN);
			pca.PCA.log(" 7. Transposing");
			expressionMatrix.transpose();
			//Centering is done before the PCA inside the PCA function, to assure the centering is done over the right dimension (PCA function also puts matrixes in right orientation)
			pca.PCA.log(" 8. Calculating column averages");
			Matrix colAverages = expressionMatrix.calcAvgCols();//NOTE: this is a vector that is based on rows already being corrected for averages.
			pca.PCA.log(" 9. Centering: Adjusting for column averages");
			expressionMatrix.adjustForAverageAllCols(colAverages);
			colAverages.write(writeFolder+ "SAMPLE_QuantNorm_columnAverages.txt");
			pca.PCA.log("10. Calculating row averages");
			Matrix rowAverages = expressionMatrix.calcAvgRows();
			//rowAverages.write(writeFolder+ "SAMPLE_QuantNorm_rowAverages.txt");
			pca.PCA.log("11. Adjusting for row averages");
			//expressionMatrix.adjustForAverageAllrows(rowAverages);
			//String expNormLogCentFile= writeFolder+"MATRIX_CenteredOverRowsAndColumns.txt";
			pca.PCA.log("12. Writing centered file in: " + expNormLogCentFile);
			expressionMatrix.write(expNormLogCentFile);
			
			//this function is slow ...
			pca.PCA.log("13. creating covariance matrix over the samples");
			covMat = pca.PCA.covariance(Paths.get(expNormLogCentFile), writeFolder+"SAMPLE_covariance.txt");
			//Use Marc Jan's script instead...
		    //java  -Xmx30g -jar javaPCoA-1.0.4-SNAPSHOT-jar-with-dependencies.jar -m covariation -r -t 20 -i /Volumes/Promise_RAID/GeneNetwork/Sipko/est_counts_nocancernocellline/MATRIX_CenteredOverRowsAndColumns.txt
			//gunzip MATRIX_CenteredOverRowsAndColumns.CovariationMatrix.txt.gz
		}  	
		else
		{
			pca.PCA.log("13.1. Reading covariance matrix");
			covMat = new MatrixStruct(writeFolder+"MATRIX_CenteredOverRowsAndColumns.CovariationMatrix.txt");
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
