package Analyses;

import java.io.IOException;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.TransformerException;

import PCA.CorrelationLargeGenes;
import PCA.CreateGeneEigenvectorFile;
import PCA.Matrix;
import PCA.MatrixStruct;
import PCA.Var;
import no.uib.cipr.matrix.NotConvergedException;

public class CorrelationVSexpressionPlots {
	//creates a file that contains in the first column the absolute average expression, in the 2nd column the average expression, 3rd column stdev in the normalized expression file
	//requires a normalized exrpession file as input
	static int nThreads = 16;
	
	public static void main(String[] args) throws IOException, NotConvergedException, InterruptedException, ParserConfigurationException, TransformerException
	{
		CreateGeneEigenvectorFile.var = new Var();
		
		//args = new String[]{"json=E:/Groningen/Data/PublicSamples/09-2016/31.07.pc1.illumina.genes.expressed.DEseqnorm/config.json"};
		
		CreateGeneEigenvectorFile.checkArgs(args);//should write one for this class itself but cba...
		if( CreateGeneEigenvectorFile.var.writeFolder == null)
			 CreateGeneEigenvectorFile.var.writeFolder = CreateGeneEigenvectorFile.var.expFile.replace(".txt", "").replace(".gz", "")+"/";
		
		CreateGeneEigenvectorFile.writeParameters();

		CreateGeneEigenvectorFile.var.filePathsExist();

		//normalize and center the data
		CreateGeneEigenvectorFile.normalize();
		
		String preCorrelationFile = CreateGeneEigenvectorFile.var.tempName;
		String type = "covariance";
		if(CreateGeneEigenvectorFile.var.correlation)	
			type = "correlation";
		String covMatFN = CreateGeneEigenvectorFile.var.writeFolder+"gene_"+type+".txt";
		Matrix preCorrelation = new Matrix(preCorrelationFile);
		CorrelationLargeGenes.correlation(covMatFN, preCorrelation, CreateGeneEigenvectorFile.var.correlation, nThreads);
		
		MatrixStruct avgExpression = new MatrixStruct(CreateGeneEigenvectorFile.var.writeFolder+ "SAMPLE_Norm_GeneAverages.txt");
		MatrixStruct avgAbsCorrelation = new MatrixStruct(covMatFN.replace(".txt", "_absoluteAverages.txt"));
		MatrixStruct stDevs = new MatrixStruct(CreateGeneEigenvectorFile.var.writeFolder + "gene_STDevs.txt");
		
		avgExpression=avgExpression.mergeColumns(avgAbsCorrelation);
		avgExpression=avgExpression.mergeColumns(stDevs);
		
		avgExpression.write(CreateGeneEigenvectorFile.var.writeFolder + "avgExpressionVSabsCorrelation.txt");
	}

}
