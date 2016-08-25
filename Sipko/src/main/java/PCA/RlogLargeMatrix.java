package PCA;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

public class RlogLargeMatrix 
{
	//This is doing the DESeq normalization, which does not use replicate information.
	//Also works for matrixes of larger then MAXint sizes
	
	//static String expressionFN = "E:/Groningen/Data/Annique/LLD_and_BIOS_Kallisto_EstimatedTranscriptCounts.ProbesWithZeroVarianceRemoved_removedBadSamples.txt.gz";
	static String expressionFN = "E:/Groningen/Data/Juha/Calculon/JuhaMerged/31.txt";//"/groups/umcg-wijmenga/tmp04/umcg-jkarjalainen/31.txt";
	static String writeFolder = "E:/Groningen/Data/Juha/Calculon/JuhaMerged/test/";//if null becomes new File(expressionFN).getParent()+"/";
	static String geoFN = null;//if null calculates geometric means based on this dataset
	static boolean writeAll = true;	//write all intemediary files too
	static double logAdd = 0.5;//previously used to multiply the results by this number, but seems pointless since it does not have any effect, so does not do anything anymore
	
	public static void main(String[] args) throws IOException 
	{
		checkArgs(args);
		if(writeFolder == null)
			writeFolder = new File(expressionFN).getParent()+"/DESeqNorm/";
		if(!new File(writeFolder).exists())
		{
			new File(writeFolder).mkdir();
		}
		if(geoFN == null)
			geoFN = writeFolder+ "geoMean.txt";
		
		Matrix expression = new Matrix(expressionFN);
		
		double start = System.nanoTime();
		rLog(writeFolder, expression, writeAll, geoFN);
		expression.logTransform(2,logAdd);//adds +0.5 before log
		double end = System.nanoTime();
		System.out.println((end-start)/1000/1000 + " ms");
		expression.write(writeFolder + new File(expressionFN).getName().replace(".txt", "").replace(".gz","")+".DESeqNorm.Log2.txt.gz");
	}

	public static void rLog(String writeFolder, Matrix expression, boolean writeAll, String writeGeoFN) throws IOException 
	{
		expression.putGenesOnRows();
		String swapFN = writeFolder + "swapFile.txt";
		expression.write(swapFN);
		JuhaPCA.PCA.log(" 6. Rlog without log");
		String correctedNotLogged =  writeFolder + new File(expressionFN).getName().replace(".txt", "").replace(".gz","") + ".DESeqNorm.txt.gz";
		rLog(expression, writeFolder, swapFN, writeGeoFN);
		
		if(writeAll)
			expression.write(correctedNotLogged);
	}
	public static void rLog(Matrix expression, String writeFolder, String fileName, String writeGeoFN) throws IOException 
	{
		rLog(expression, writeFolder, fileName, null, writeGeoFN);
	}
	public static void rLog(Matrix expression, String writeFolder, String fileName, MatrixStruct geoMean, String writeGeoFN) throws IOException 
	{
		if(geoMean == null)
		{
			geoMean = getGeoMeans(expression, writeFolder, writeGeoFN);
			//need to read the matrix again here...
			expression.readFile(fileName);
		}

		Matrix geoMeanMat = new Matrix(geoMean);
		geoMeanMat.keepRows(expression);//keep only the rows that are also in the geomean file
		
		expression.write(fileName.replace(".txt", "geoRowsOnly.txt"));
		MatrixStruct denominators = new MatrixStruct(expression.cols(),1);
		denominators.setRowHeaders(expression.getColHeaders());
		
		//determine the denominator per sample
		for(int c = 0; c < expression.cols(); c++)
		{
			//get the correct denominator
			double[] column = new double[expression.rows()];
			for(int r = 0; r < column.length; r++)
			{
				column[r] = (expression.matrix.get(r,c)) / geoMean.matrix.get(r,0);
			}
			Arrays.sort(column);
			
			org.apache.commons.math3.stat.descriptive.rank.Median med = new org.apache.commons.math3.stat.descriptive.rank.Median();
			denominators.matrix.set(c,0,med.evaluate(column));
		}
		if(writeFolder != null)
			denominators.write(writeFolder + "Denominators.txt");
		
		expression.readFile(fileName);
		//calculate the normalized readcounts
		for(int c = 0; c < expression.cols(); c++)
		{
			for(int r = 0; r < expression.rows(); r++)
			{
				expression.matrix.set(r,c,expression.matrix.get(r,c)/denominators.matrix.get(c,0));	
			}
		}
	}
	
	private static MatrixStruct getGeoMeans(Matrix expression ,String writeFolder, String writeGeoFN) throws IOException {
		expression.logTransform(10,0);
		MatrixStruct geoMean = new MatrixStruct(expression.rows(),1);
		geoMean.setRowHeaders(expression.getRowHeaders());
		MatrixStruct keepGenes = new MatrixStruct(expression.rows(),1);
		for(int r =0; r < expression.rows(); r++)
		{
			double gM = expression.sumRow(r)/expression.cols();
			if(Double.isFinite(gM))
			{
				geoMean.matrix.set(r,0,gM);
				geoMean.matrix.set(r,0,Math.pow(10,geoMean.matrix.get(r,0)));
				keepGenes.setRowHeader(r, expression.getRowHeaders()[r]);
			}
		}
		
		if(writeFolder != null)
		{
			geoMean.keepRows(keepGenes);
			System.out.println("geofilename=" + writeGeoFN);
			geoMean.write(writeGeoFN);
		}
		return geoMean;
	}
	
	static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko") && args.length < 1)
			return;
		if(args.length < 1)
		{
			System.out.println("Script requires the following argumetns:\n"
					+ "1. filename=<expressionFN.txt> - Expression file with genes on rows samples on columns\n"
					+ "2. writeFolder=<writeFolderFN.txt> - Folder where the files will be written (default=parentFolder(input.txt))\n"
					+ "3. geoFN=<geoFn.txt> - Optional file with geometric mean per gene to use (default=null)\n");
			System.exit(1);
		}
		
		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase())
			{
				case "filename":
					expressionFN =value;
					break;
				case "writefolder":
					writeFolder = value;
					break;
				case "geofn":
					geoFN = value;
					break;
				case "writeall":
					writeAll = Boolean.parseBoolean(value);
					break;
				case "totalreadcount":
					logAdd = Double.parseDouble(value);
				break;
				default:
					System.out.println("Incorrect argument supplied:\n"+ args[a] +"\nexiting");
					System.exit(1);
			}
		}
	}
}
