package PCA;

import java.io.IOException;

import pca.MatrixStruct;

public class HighestExpressed 
{
	public static void main(String[] args) throws IOException 
	{
		String expressionFN = "";
		String writeFolder = "";
		boolean writeAll = true;
		boolean tpm = false;
		double correctTotalReadCount = 0;
		double highestExpressed =1;
		boolean STdevCutoff = false;
		writeAll = true;
		
		
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
				case "stdevcutoff":
					STdevCutoff = Boolean.parseBoolean(value);
					break;
				case "tpm":
					tpm = Boolean.parseBoolean(value);
					break;
				case "correcttotalreadcount":
					correctTotalReadCount = Double.parseDouble(value);
					break;
				case "highestexpressed":
					highestExpressed = Double.parseDouble(value);
					break;
				default:
					checkArgs(args);
					System.out.println("Incorrect argument supplied; exiting");
					System.exit(1);
			}
		}
		
		MatrixStruct expressionStruct = new MatrixStruct(expressionFN);
		highestExpressed(expressionStruct, tpm, correctTotalReadCount, writeFolder, highestExpressed, STdevCutoff, writeAll);
	}

	public static void highestExpressed(MatrixStruct expressionStruct, boolean tpm, double correctTotalReadCount,
			String writeFolder, double highestExpressed, boolean STdevCutoff, boolean writeAll) throws IOException 
	{
		if(!tpm && correctTotalReadCount <1)
		{
			pca.PCA.log(" Quantile normalization before taking averageCutoff");
			expressionStruct.expressionToRank(expressionStruct.quantileNormVector(), 0);
		}	
		String averagesFN = writeFolder+"AveragesAllGenes.txt";
		expressionStruct.getAveragesPerRow(expressionStruct)
						.write(averagesFN);
		pca.PCA.log("  . Removing the " + ((1.0-highestExpressed)*100) + " percent lowest expressed genes" );
		
		if(STdevCutoff)
		{
			MatrixStruct stDevs = expressionStruct.stDevRows();
			String stdevFN = writeFolder + "gene_STDevsForCutoff.txt";
			stDevs.write(stdevFN);
			PCcorrection.keepTopPercentage(expressionStruct,stdevFN, highestExpressed, averagesFN.replace(".txt", "top_" + highestExpressed+ ".txt"), true, writeAll);
		}
		else
			PCcorrection.keepTopPercentage(expressionStruct,averagesFN, highestExpressed, averagesFN.replace(".txt", "top_" + highestExpressed+ ".txt"), false, writeAll);
	}
	static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
			return;
		System.out.println("Wrong arguments");
		System.exit(1);
	}
}
