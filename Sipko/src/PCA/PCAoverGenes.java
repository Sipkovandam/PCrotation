package PCA;

import java.io.IOException;
import java.nio.file.Paths;

import no.uib.cipr.matrix.NotConvergedException;
import pca.MatrixStruct;
import umcg.genetica.math.stats.concurrent.ConcurrentCorrelation;
import umcg.genetica.math.stats.concurrent.ConcurrentCovariation;

public class PCAoverGenes 
{
	public static void main(String[] args) throws IOException, NotConvergedException
	{
		String expFile = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/" + "TESTexpression.txt";
		String chromLocationsFile = "E:/Groningen/Data/GenePositionInfo_23X_24Y_25MT_26rest.txt";
		//String expFile = "E:/Groningen/Data/PublicSamples/ComputerTest/Rtest/" + "rSample.txt";
		System.out.println(System.getProperty("user.dir"));
		CreateGeneEigenvectorFile.checkArgs(args);
		
		boolean writeAll = true;
		
		if(!System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko") || args.length !=0)
		{
			expFile = args[0];
			if(args.length>1)
				chromLocationsFile = args[1];
			else
				chromLocationsFile = null;
			if(args.length>2 && args[2].contains("nowrite"))
				writeAll = false;
		}
		
		String writeFolder = expFile.replace(".txt", "_OverGenes/");
		CreateGeneEigenvectorFile.makeFolder(writeFolder);
		
		MatrixStruct covMat = null;
		String expNormLogCentFile= writeFolder+"MATRIX_Centered.txt";
		MatrixStruct expressionMatrixStruct = new MatrixStruct(expFile);
		
		if(chromLocationsFile != null)
		{
			pca.PCA.log(" 0. Sorting IDs based on chromosome locations");
			pca.PCA.log(" 1. Transposing");
			expressionMatrixStruct.transpose();
			MatrixStruct chromLocs = new MatrixStruct(chromLocationsFile);
			expressionMatrixStruct = SortChromosome.sort(expressionMatrixStruct, chromLocs);
		}
		else
		{
			pca.PCA.log(" 1. Transposing");
			expressionMatrixStruct.transpose();
		}
		expressionMatrixStruct.write(writeFolder+"startMatrix.txt");
		pca.PCA.log(" 1.1 Adding random values <5 to values <5");
		double randomAddValue = 5;
		expressionMatrixStruct.addRandomValues(randomAddValue);
		expressionMatrixStruct.write(writeFolder+"randAddedToBelow_" +randomAddValue + ".txt");
		Matrix expressionMatrix = null;
		expressionMatrix = new Matrix(expressionMatrixStruct);
		expressionMatrixStruct = null;
		
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
		pca.PCA.log(" 9. Calculating column averages");
		Matrix colAverages = expressionMatrix.calcAvgCols();//NOTE: this is a vector that is based on rows already being corrected for averages.
		colAverages.write(writeFolder+ "SAMPLE_QuantNorm_columnAverages.txt");
		pca.PCA.log("10. Centering: Adjusting for column averages");
		expressionMatrix.adjustForAverageAllCols(colAverages);
		expressionMatrix.write(quantFN.replace(".txt", "_adjustedForGeneAverages.txt"));
		pca.PCA.log("10. Calculating row averages");
		Matrix rowAverages = expressionMatrix.calcAvgRows();
		rowAverages.write(writeFolder+ "SAMPLE_QuantNorm_rowAverages.txt");
		pca.PCA.log("11. Adjusting for row averages");
		expressionMatrix.adjustForAverageAllrows(rowAverages);
		//String expNormLogCentFile= writeFolder+"MATRIX_CenteredOverRowsAndColumns.txt";
		pca.PCA.log("11. Writing centered file in: " + expNormLogCentFile);
		expressionMatrix.write(expNormLogCentFile);
		/**6. Calculate PCscores of the genes of public/healthy expression data**/
		
		MatrixStruct expressionStruct = new MatrixStruct(expressionMatrix.rowNames, expressionMatrix.colNames, expressionMatrix.values);
		pca.PCA.log("12. creating covariance matrix over the GENES");
		ConcurrentCorrelation calculator = new ConcurrentCorrelation(20);//correlation instead of covariation
		//lets do it over the genes! ;)
		expressionStruct.transpose();
		double[][] inMat = expressionStruct.getSquareMatrix();
		double[][] covMatrix = calculator.pairwiseCorrelation(inMat);;//correlation instead of covariation
		covMat = new MatrixStruct(expressionStruct.getRowHeaders(), expressionStruct.getRowHeaders(), covMatrix);
		covMat.write(writeFolder+"GENE_covariance.txt");
		inMat = null; covMatrix = null; System.gc();System.gc();
		
		pca.PCA.log("14. calculating eigenvalues over the GENES");
		MatrixStruct[] evds = pca.PCA.evd(covMat, Paths.get(writeFolder+"GENE"));
		MatrixStruct eigenVectors = evds[0];
		MatrixStruct PCeigenvalues = evds[1];
	}
}
