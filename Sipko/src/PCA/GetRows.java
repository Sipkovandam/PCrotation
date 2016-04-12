package PCA;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;

public class GetRows 
{

	public static void main(String args[]) throws IOException
	{
		String fileName = "E:/Groningen/Data/PublicSamples/Test13/TPM_9900Samples/ranks.txt";
		//String fileName2 = "E:/Groningen/Data/PublicSamples/Test9/PublicSamplesWithoutDownSyndrome.txt";
		String fileName2 = "ENSG00000268903,ENSG00000269981,ENSG00000225630";
		String writeName = "E:/Groningen/Data/PublicSamples/Test13/TPM_9900Samples/ranks_small.txt";
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
						fileName2 = value;
						break;
					case "rowstoget":
						fileName2 = value;
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
		
		if(!fileName2.contains(",") && fileName2.contains(".txt"))
			file2 = new Matrix(fileName2);
		else{
			String[] rowNames = fileName2.split(",");
			file2 = new Matrix(rowNames.length, 1);
			file2.rowNames = rowNames;
			file2.colNames[0] = "-";
		}
		
		BufferedReader reader = new BufferedReader(new FileReader(new File(fileName)));
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(writeName)));
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
					+ "It uses the following 2 arguments:\n"
					+ "1. fileName=<fileName> -  File to retrieve rows from\n"
					+ "2.1 rowstoget=<fileName> -  File with the rows to keep in the first column(header row is ignored)\n"
					+ "2.2 rowstoget=<gene1,gene2,gene3> - Alternatively you can use a comma separated list of rowNames you wish to retrive (SRR001,SR002,...)\n"
					+ "3. writeName=<fileName> - Name of the file to write to \n"
					+ "Make sure each file has at least 2 columns and rows (just at a bunch of 0's in the 2nd column if you must) \n");
			System.exit(1);
		}
	}
}
