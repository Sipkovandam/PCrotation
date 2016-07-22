package PCA;

import java.io.IOException;

import eqtlmappingpipeline.normalization.Normalizer;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

public class GCcontent 
{
	static String expressionFN = "";
	static String gcContentFN = "";
	static String writeFN = null;

	public static void main (String[] args) throws Exception
	{
		System.out.println(gcContentFN);
		checkArgs(args);
		System.out.println(gcContentFN);
		MatrixStruct geneGC = new MatrixStruct(gcContentFN);
		
		MatrixStruct expression = new MatrixStruct(expressionFN);
		
		calculateAndCorrect(expression, geneGC, "GCperSample.txt", "GC_corrected.txt");	
	}
	static MatrixStruct calculateGCperSample(MatrixStruct expression, MatrixStruct geneGC, String writeFN) {
		expression.putGenesOnRows();
		MatrixStruct sampleGC = new MatrixStruct(expression.getColHeaders().length,1);
		try
		{
			sampleGC.setRowHeaders(expression.getColHeaders());
			sampleGC.setColHeaders(new String[]{"GC content"});
			
			for(int c = 0; c < expression.cols(); c++)
			{
				double totalGC = 0;
				double totalReads = 0;
				for(int r = 0; r < expression.rows(); r++)
				{
					String gene = expression.getRowHeaders()[r];
					if(geneGC.rowHash.get(gene)== null)
						continue;
					double gcPercentage =  geneGC.matrix.get(geneGC.rowHash.get(gene), 0)/100;
					//System.out.println(gcPercentage);
					totalGC += expression.matrix.get(r, c)*gcPercentage;
					//System.out.println(totalGC);
					totalReads += expression.matrix.get(r, c);
				}
				//System.out.println("GCcontent =" + totalGC+"/"+totalReads + "=" + totalGC/totalReads);
				sampleGC.matrix.set(c, 0, totalGC/totalReads);
			}
			sampleGC.write(writeFN);
			System.out.println("Done, file written to: " + writeFN);
		}catch (Exception e){e.printStackTrace();}
		return sampleGC;
	}
	public static MatrixStruct calculateAndCorrect(MatrixStruct expressionStruct, MatrixStruct geneGC, String gCperSampleWriteFN, String gcCorrectedWriteFN)
	{
		GCcontent.calculateGCperSample(expressionStruct, geneGC,gCperSampleWriteFN);
		
		Normalizer normalizer = new Normalizer();
		System.out.println("rows = " + expressionStruct.rows() + " cols = " + expressionStruct.cols() );
		System.out.println(expressionStruct.getRowHeaders()[3]);
		DoubleMatrixDataset<String,String> data = new DoubleMatrixDataset<String,String>(expressionStruct.getMatrix(), expressionStruct.getRowHeadersAsList(),expressionStruct.getColHeadersAsList());
		try {
			normalizer.adjustCovariates(data, gcCorrectedWriteFN, gCperSampleWriteFN, false, 1E-10);
		} catch (IOException e) {e.printStackTrace();}
		//this is needed because expressionStruct.getMatrix (above) is a deepcopy
		return new MatrixStruct(expressionStruct.getRowHeaders(), expressionStruct.getColHeaders(), data.getRawData());
	}
	static void checkArgs(String[] args) 
	{
		if(args.length < 1)
		{
			System.out.println("Wrong arguments, needs:\n"
					+ "1. expressionFN=<expressionFile.txt> - filename of the file containing the expression per gene per sample \n"
					+ "2. gcContentFN=<gcContentFile.txt> - filename of the file containing the GC percentage (e.g. 40.12) per gene (1st column gene name, 2nd column GC content (without % sign)\n");
			System.exit(1);
		}
		
		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase())
			{
				case "expressionfn":
					expressionFN =value;
					break;
				case "gccontentfn":
					gcContentFN =value;
					break;
				default:
					System.out.println("Incorrect argument supplied:\n"+ args[a] +"\nexiting");
					System.exit(1);
			}
		}
		if(writeFN == null)
			writeFN = expressionFN.replace(".txt", "_GCcontent.txt");
	}
}
