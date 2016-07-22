package PCA;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;

import umcg.genetica.math.stats.concurrent.ConcurrentCovariation;

public class KeepTopVariance 
{
	public static void main(String[] args) throws IOException 
	{
		String expressionFN = "";
		String writeFolder = "";
		boolean writeAll = true;	
		boolean correl = false;
		int ignoreLowestValues = -1;
		double topVariance = 1;

		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase()){
				case "filename":
					expressionFN =value;
					break;
				case "writefolder":
					writeFolder = value;
					break;
				case "writeall":
					writeAll = Boolean.parseBoolean(value);
					break;
				case "correl":
					correl = Boolean.parseBoolean(value);
					break;
				case "ignorelowestvalues":
					ignoreLowestValues = Integer.parseInt(value);
					break;
				case "topvariance":
					topVariance = Double.parseDouble(value);
					break;
				default:
					checkArgs(args);
					System.out.println("Incorrect argument supplied; exiting");
					System.exit(1);
			}
		}
		
		MatrixStruct expressionStruct = new MatrixStruct(expressionFN);
		keepTopVariance(expressionStruct, correl, ignoreLowestValues, writeFolder, topVariance);
	}

	public static void keepTopVariance(MatrixStruct expressionStruct, boolean correl, int ignoreLowestValues, String writeFolder, double topVariance) throws IOException 
	{
		//create a covariance matrix over the genes
		//calculate the average absolute covariance per gene
		pca.PCA.log("21.1 creating covariance/correlation matrix OVER THE GENES");

		//genes are on x-Axis after this:
		expressionStruct.transpose();
		
		ConcurrentCovariation calculatorGenes = new ConcurrentCovariation(20);
		double[][] inMat = expressionStruct.getMatrix();		
		
		String type = "";
		if(correl)
			type = "correlation";
		else
			type = "covariance";
			
		double[] cutoffs = getCutoffs(ignoreLowestValues, expressionStruct); 
		if(cutoffs!= null)
			for(double c : cutoffs)
				System.out.println(c);
		calculatorGenes.pairwiseCovariation(inMat,true, writeFolder+"gene_"+type+".txt", expressionStruct.getRowHeaders(), correl, cutoffs);//last argument, if false = covariance; true = correlation.
		inMat = null; System.gc(); System.gc();
		System.out.println("topvar = " + topVariance + " expressionStruct.rows();" + expressionStruct.rows());
		
		PCcorrection.keepTopPercentage(expressionStruct,writeFolder+"gene_"+type+"_Average.txt", topVariance, writeFolder+"gene_"+type+"_remainingGenes.txt", true, true);
		System.out.println("topvar = " + topVariance + " expressionStruct.rows();" + expressionStruct.rows());
		expressionStruct.transpose();
	}
	public static double[] getCutoffs(int removeLowestValues, MatrixStruct genesOnRows) {
		if(removeLowestValues < 1)
			return null;
		double[] cutoffs = new double[genesOnRows.cols()];
		for(int c = 0; c < genesOnRows.cols(); c++){
			cutoffs[c] = Double.POSITIVE_INFINITY;
			Hashtable<Double, Integer> values = new Hashtable<Double,Integer>();
			ArrayList<Double> vals = new ArrayList<Double>();
			for(int r = 0; r < genesOnRows.rows(); r++){
				if(!values.containsKey(genesOnRows.matrix.get(r, c))){
					vals.add(genesOnRows.matrix.get(r, c));
					values.put(genesOnRows.matrix.get(r, c), r);
				}
			}
			Collections.sort(vals);;
			int get = removeLowestValues-1;
			if(removeLowestValues > vals.size())
				get = vals.size()-1;
			cutoffs[c] = vals.get(get);
		}
		
		return cutoffs;
	}
	static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
			return;
		System.out.println("Wrong arguments");
		System.exit(1);
	}
}
