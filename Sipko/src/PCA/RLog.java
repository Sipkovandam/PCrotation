package PCA;

import java.io.File;
import java.io.IOException;

import pca.MatrixStruct;

public class RLog 
{
	static String expressionFN = "E:/Groningen/Data/Annique/LLD_and_BIOS_Kallisto_EstimatedTranscriptCounts.ProbesWithZeroVarianceRemoved_removedBadSamples.txt.gz";
	static String writeFolder = null;//if null becomes new File(expressionFN).getParent()+"/";
	static String geoFN = null;//if null calculates geometric means based on this dataset
	static boolean writeAll = true;	//write all intemediary files too
	static double rLog = 1;//previously used to multiply the results by this number, but seems pointless since it does not have any effect, so does not do anything anymore
	
	public static void main(String[] args) throws IOException 
	{
		checkArgs(args);
		if(writeFolder == null)
			writeFolder = new File(expressionFN).getParent()+"/";
		if(geoFN == null)
			geoFN = writeFolder+ "geoMean.txt";
		

		MatrixStruct expressionStruct = new MatrixStruct(expressionFN);
		rLog(writeFolder, expressionStruct, rLog, writeAll, geoFN);
		expressionStruct.log2Transform(1);
		expressionStruct.write(expressionFN.replace(".txt", "").replace(".gz",  "") + ".deSeqNorm.Log2.txt.gz");
	}

	public static void rLog(String writeFolder, MatrixStruct expressionStruct, double rLog, boolean writeAll, String writeGeoFN) throws IOException 
	{
		String swapFN = writeFolder + "swapFile.txt";
		expressionStruct.write(swapFN);
		pca.PCA.log(" 6. Rlog without log");
		String correctedNotLogged =  writeFolder+ ".deSeqNorm.txt.gz";
		expressionStruct.rLog(rLog, writeFolder, swapFN, writeGeoFN);
		if(writeAll)
			expressionStruct.write(correctedNotLogged);
	}
	static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
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
					rLog = Double.parseDouble(value);
				break;
				default:
					System.out.println("Incorrect argument supplied:\n"+ args[a] +"\nexiting");
					System.exit(1);
			}
		}
	}
}
