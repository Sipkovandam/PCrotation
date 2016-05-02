package PCA;

import java.io.IOException;

import pca.MatrixStruct;
import umcg.genetica.math.stats.concurrent.ConcurrentCovariation;

public class CorrelationMatrixLarge {

	public static void main(String[] args) throws IOException 
	{
		String expressionFN = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/TESTexpression.txt";
		String writeFile = null;
		boolean correlation = false;
		
		if(args.length==0) checkArgs(args);
		
		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase()){
				case "filename":
					expressionFN =value;
					break;
				case "writefile":
					writeFile = value;
					break;
				case "correlation":
					correlation = Boolean.parseBoolean(value);
					break;
				default:
					checkArgs(args);
					System.out.println("Incorrect argument supplied; exiting");
					System.exit(1);
			}
		}
		if(writeFile == null)
		{
			if(correlation)
				writeFile = expressionFN.replace(".txt", "_correlation.txt");
			else//if covariance
				writeFile = expressionFN.replace(".txt", "_covariance.txt");
		}	
		Matrix expression = new Matrix(expressionFN);
		
		correlation(writeFile, expression,correlation);
		System.out.println("Done, file written to: "+ writeFile);
	}

	public static void correlation(String writeFile, Matrix expression, boolean correlation) throws IOException 
	{
		expression.putGenesOnRows();
		
		System.out.println("Centering data over rows");
		Matrix rowAverages = expression.calcAvgRows();
		rowAverages.write(writeFile.replace(".txt", "_inputRowAverages.txt"));
		expression.adjustForAverageAllrows(rowAverages);
		
		ConcurrentCovariation calculatorGenes = new ConcurrentCovariation(20);
		double[][] matrix = calculatorGenes.pairwiseCovariation(expression.values,false, writeFile, expression.rowNames,correlation, null);//last argument, if false = covariance; true = correlation.
		Matrix covMat = new Matrix(expression.rowNames, expression.rowNames, matrix);
		covMat.write(writeFile);
	}
	static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
			return;
		System.out.println("This script uses the following arguments:\n"
				+ "1. fileName=<fileName> \n"
				+ "2. writeFile=<writeFileName> default(covariance.txt/correlation.txt)\n"
				+ "3. correlation=<false/true> default(true)"
				+ "example:\n"
				+ "java -jar -Xmx50G CorrelationLarge.jar fileName=expression.txt correlation=false");
		System.exit(1);
	}
}
