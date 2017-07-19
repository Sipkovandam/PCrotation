package PCA;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import MatrixScripts.MyMatrix;
import Tools.Script;

public class DeSeqNorm
{
	//This is doing the DESeq normalization, which does not use replicate information.
	
	//static String expressionFN = "E:/Groningen/Data/Annique/LLD_and_BIOS_Kallisto_EstimatedTranscriptCounts.ProbesWithZeroVarianceRemoved_removedBadSamples.txt.gz";
	static String expressionFN = "E:/Groningen/Scripts/Tests/Rlog.java/Samples.txt";
	
	static String writeFolder = null;//if null becomes new File(expressionFN).getParent()+"/";
	static String geoFN = null;//if null calculates geometric means based on this dataset
	static boolean writeAll = true;	//write all intemediary files too
	static double logAdd = 0.5;//previously used to multiply the results by this number, but seems pointless since it does not have any effect, so does not do anything anymore
	static boolean log = false;
	static boolean genes = true;
	static boolean roundValues = true;//rounds expression values to whole counts
	
	public static void main(String[] args) throws IOException 
	{
		checkArgs(args);
		if(writeFolder == null)
		{
			writeFolder = new File(expressionFN).getParent()+"/DESeqNorm/";
			
		}
		if(!new File(writeFolder).exists())
			new File(writeFolder).mkdirs();
		
		MyMatrix expressionStruct = new MyMatrix(expressionFN);
		if(genes)
			expressionStruct.putGenesOnRows();
		if(roundValues)
			expressionStruct.roundValues();
		double start = System.nanoTime();
		rLog(writeFolder, expressionStruct, writeAll, geoFN);
		double end = System.nanoTime();
		System.out.println((end-start)/1000/1000/1000 + " sec");
		if(log)
		{
			expressionStruct.log2Transform(logAdd);
			//write the results
			expressionStruct.write(writeFolder + new File(expressionFN).getName().replace(".txt", "").replace(".gz","")+".DESeqNorm.Log2.txt.gz");
		}
	}

	public static void rLog(String writeFolder, MyMatrix expressionStruct, boolean writeAll, String geoFN) throws IOException 
	{
		String swapFN = writeFolder + "swapFile.txt";
		expressionStruct.write(swapFN);
		JuhaPCA.PCA.log(" 6. Deseq normalization without log");
		String correctedNotLogged =  writeFolder + new File(expressionFN).getName().replace(".txt", "").replace(".gz","") + ".DESeqNorm.txt.gz";
		if(geoFN!= null && new File(geoFN).exists())
			rLog(expressionStruct, writeFolder, swapFN, new MyMatrix(geoFN), null);
		else
			rLog(expressionStruct, writeFolder, swapFN, null, geoFN);
		if(writeAll)
			expressionStruct.write(correctedNotLogged);
	}
	
	public static void rLog(MyMatrix expressionStruct, String writeFolder, String fileName, String writeGeoFN) throws IOException 
	{
		rLog(expressionStruct, writeFolder, fileName, null, writeGeoFN);
	}
	
	private static MyMatrix getGeoMeans(MyMatrix expressionStruct, String writeFolder, String writeGeoFN) throws IOException {
		//log transform data so values do not become to large (and you save a lot of speed (you can add instead of multiply))
		expressionStruct.logTransform(10,0);
		//create the matrix where the results are stored
		MyMatrix geoMean = new MyMatrix(expressionStruct.rows(),1);
		geoMean.setRowHeaders(expressionStruct.getRowHeaders());
		//Creates a matrix where you can put the genes you want to keep (would normally use a hashtable, but this works with MatrixStruct.keeprows called later)
		MyMatrix keepGenes = new MyMatrix(expressionStruct.rows(),1);
		for(int r =0; r < expressionStruct.rows(); r++)
		{
			//Calculate the geometric mean (based on the logged values, which is why it is summing the row instead of multiplying. Also why it is dividing instead of ^(1/x))
			double gM = expressionStruct.sumRow(r)/expressionStruct.cols();
			if(Double.isFinite(gM))
			{
				geoMean.matrix.set(r,0,gM);
				geoMean.matrix.set(r,0,Math.pow(10,geoMean.matrix.get(r,0)));//raise 10^X to unlog the resulting values again
				keepGenes.setRowHeader(r, expressionStruct.getRowHeaders()[r]);
			}
		}
		
		geoMean.keepRows(keepGenes); 
		if(writeGeoFN==null)
			writeGeoFN = writeFolder+"geoMean.txt";
		System.out.println("geofilename=" + writeGeoFN);
		geoMean.write(writeGeoFN);
		
		return geoMean;
	}
	
	public static void rLog(MyMatrix expressionStruct, String writeFolder, String fileName, MyMatrix geoMean, String writeGeoFN) throws IOException 
	{
		if(geoMean == null)
		{
			geoMean = getGeoMeans(expressionStruct, writeFolder, writeGeoFN);
			//need to read the matrix again here because we logged it before...
			expressionStruct.readFile(fileName);
		}
		MyMatrix denominators = getDenominators(expressionStruct, fileName, geoMean, writeFolder);
		
		expressionStruct.readFile(fileName);//read the file again since we kept only the rows that have no 0 values.
		
		//calculate the normalized readcounts
		for(int c = 0; c < expressionStruct.cols(); c++)
			for(int r = 0; r < expressionStruct.rows(); r++)
				expressionStruct.matrix.set(r,c,expressionStruct.matrix.get(r,c)/denominators.matrix.get(c,0));	
	}
	
	private static MyMatrix getDenominators(MyMatrix expressionStruct, String fileName, MyMatrix geoMean, String writeFolder) throws IOException {
		geoMean.keepRows(expressionStruct);//keep only the rows that are also in the geomean file
		expressionStruct.write(fileName.replace(".txt", "geoRowsOnly.txt"));
		MyMatrix denominators = new MyMatrix(expressionStruct.cols(),1);
		denominators.setRowHeaders(expressionStruct.getColHeaders());
		
		//determine the denominator per sample
		for(int c = 0; c < expressionStruct.cols(); c++)
		{
			//get the correct denominator
			double[] column = new double[expressionStruct.rows()];
			for(int r = 0; r < column.length; r++)
			{
				column[r] = (expressionStruct.matrix.get(r,c)) / geoMean.matrix.get(r,0);
			}
			Arrays.sort(column);
			
			org.apache.commons.math3.stat.descriptive.rank.Median med = new org.apache.commons.math3.stat.descriptive.rank.Median();
			denominators.matrix.set(c,0,med.evaluate(column));
		}
		if(writeFolder != null)
			denominators.write(writeFolder + "Denominators.txt");
		return denominators;
	}

	static void printUsage()
	{
		System.out.println("Script requires the following arguments:\n"
				+ "1. filename=<expressionFN.txt> - Expression file with genes on rows samples on columns\n"
				+ "2. writeFolder=</writeFolderFN/> - Folder where the files will be written (default=parentFolder(input.txt))\n"
				+ "3. geoFN=<geoFn.txt> - Optional file with geometric mean per gene to use,\n"
				+ "   if not supplied one is created based on current dataset (default=null)\n"
				+ "4. log=<false> - If true, also a logged file will be created (default=false)\n"
				+ "5. logAdd=<0.5> - Value that is added before the log (default=0.5)\n"
				+ "Example:  java -jar Rlog.jar filename=/path/expression.txt writeFolder=/path/DESeqNorm/ log=true");
	}
	
	static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
			return;
		if(args.length < 1)
		{
			printUsage();
			System.exit(1);
		}
		
		for(int a = 0; a < args.length; a++)
		{
			String split[] = args[a].split("=");
			if(split.length<2)
			{
				printUsage();
				System.out.println("Incorrect argument supplied:\n"+ args[a] +"\nexiting");
				System.exit(1);
			}
				
			String arg = split[0];
			String value = split[1];
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
				case "logadd":
					logAdd = Double.parseDouble(value);
				break;
				case "log":
					log = Boolean.parseBoolean(value);
				break;
				case "roundvalues":
					roundValues = Boolean.parseBoolean(value);
				break;
				default:
					printUsage();
					System.out.println("Incorrect argument supplied:\n"+ args[a] +"\nexiting");
					System.exit(1);
			}
		}
	}
}
