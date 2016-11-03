package PCA;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;

import Tools.FileUtils;

public class GetRows 
{
	//gets the rows defined in the getGenesFN

	public static void main(String args[]) throws IOException
	{
		String fileName = "E:/Groningen/Data/Juha/Genes31995/Healthy/PCA/31.07.pc1.illumina.genes.expressed_DownSamples/PC_1-229_DevidedBySTdevs.txt";
		String getGenesFN = "E:/Groningen/Data/Juha/Genes31995/Old/31.07.pc1.illumina.genes.expressed.DEseqnorm_notRounded/18DownSyndrome26Normal2Cancer_counts_transposed/PC_1-300_DevidedBySTdevsTop12000Expressed.txt";
		//String fileName2 = "ENSG00000268903,ENSG00000269981,ENSG00000225630";
		String writeName = "E:/Groningen/Data/Juha/Genes31995/Healthy/PCA/31.07.pc1.illumina.genes.expressed_DownSamples/PC_1-229_DevidedBySTdevs_12000highest.txt";
		String replace = null;
		
		checkArgs(args);
		if(!System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
		{
			for(int a = 0; a < args.length; a++)
			{
				String arg = args[a].split("=")[0];
				String value = args[a].split("=")[1];
				switch (arg.toLowerCase()){
					case "filename":
						fileName = value;
						break;
					case "getgenes":
						getGenesFN = value;
						break;
					case "rowstoget":
						getGenesFN = value;
						break;
					case "writename":
						writeName = value;
						break;
					case "remove":
						replace = value;
						break;
					default:
						checkArgs(args);
						System.out.println("Incorrect argument supplied; exiting");
						System.exit(1);
				}
			}
		}
		
		
		
		Matrix file2 = null;
		
		if(!getGenesFN.contains(",") && getGenesFN.contains(".txt"))
			file2 = new Matrix(getGenesFN);
		else{
			String[] rowNames = getGenesFN.split(",");
			file2 = new Matrix(rowNames.length, 1);
			file2.rowNames = rowNames;
			file2.colNames[0] = "-";
		}
		
		BufferedReader reader = FileUtils.createReader(fileName);
		BufferedWriter writer = FileUtils.createWriter(writeName);
		String line = reader.readLine();
		writer.write(line+"\n");//write the first line by default (headers)
		Hashtable<String, Integer> toGet = file2.namesToHash(file2.rowNames);
		String[] results = new String[toGet.size()];
		while((line = reader.readLine()) != null)
		{	
			if(replace != null)
				line = line.replace(replace, "");
			String rowName = line.split("\t")[0];
			if(toGet.containsKey(rowName))
			{
				results[toGet.get(rowName)] = line;
				//writer.write(line+"\n");
			}
		}
		for(int r = 0; r < results.length; r++)
		{
			if(results[r] == null)
				continue;
			//System.out.println("r= "+(results[r]));
			writer.write(results[r]+"\n");
		}
		
		writer.close();
		reader.close();
		System.out.println("File written to:" + writeName);
	}
	public static void checkArgs(String[] args)
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
			return;
		if(args.length < 3)
		{
			System.out.println("Arguments supplied =" + args.length);
			System.out.println("This script retrieves rows for which rowNames are present in 1 file "
					+ "from another matrix file (automagically keeps header row).\n"
					+ "It uses the following 3 arguments:\n"
					+ "1. fileName=<fileName> -  File to retrieve rows from\n"
					+ "2.1 rowstoget=<fileName> -  File with the rows to keep in the first column(header row is ignored)\n"
					+ "2.2 rowstoget=<gene1,gene2,gene3> - Alternatively you can use a comma separated list of rowNames you wish to retrive (SRR001,SR002,...)\n"
					+ "3. writeName=<fileName> - Name of the file to write to \n");
			System.exit(1);
		}
	}
}
