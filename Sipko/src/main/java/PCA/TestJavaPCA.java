package PCA;

import java.io.IOException;
import java.nio.file.Paths;

import JuhaPCA.PCA;
import MatrixScripts.MatrixStruct;
import no.uib.cipr.matrix.NotConvergedException;

public class TestJavaPCA 
{
	static String writeFolder = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/TestJavaPCA/";
	static MatrixStruct covMat = new MatrixStruct("E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/TESTexpression/gene_correlation.txt");
	static MatrixStruct expressionStruct = new MatrixStruct("E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/TESTexpression/MATRIX_Centered.txt");
	
	public static void main (String args[]) throws IOException, NotConvergedException
	{
		pca();
		System.out.println("done");
	}
	public static void pca() throws IOException, NotConvergedException
	{
		MatrixStruct[] evds = JuhaPCA.PCA.evd(covMat, Paths.get(writeFolder+"gene"));
		MatrixStruct eigenVectors = evds[0];
		MatrixStruct PCeigenvalues = evds[1];
		eigenVectors.write(writeFolder+"Eigenvectors.txt");
		
//		pca.PCA.log("23. calculating PCscores over the genes");
//		String saveNamePCscoresGene = writeFolder + "GENE_PC.txt";
//		MatrixStruct[] PCscoresGenesAndAverages = PCA.scores(eigenVectors,expressionStruct, saveNamePCscoresGene);
//		MatrixStruct PCscoresGenes = PCscoresGenesAndAverages[0];
//		PCscoresGenesAndAverages=null;
	}
}
