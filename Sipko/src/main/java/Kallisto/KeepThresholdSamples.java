package Kallisto;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Set;

import MatrixScripts.MyMatrix;
import Tools.FileUtils;

public class KeepThresholdSamples 
{
	static String countsFN = "E:/Groningen/Test/JSON/ServerTest/Kallisto/Resultscounts.txt.gz";
	static String mappingPercentageFN = "E:/Groningen/Test/JSON/ServerTest/Kallisto/mappingPerSample.txt"; 
	static String writeFN = null; 
	static double threshold = 0;
	
	public static void main(String[] args) throws FileNotFoundException, IOException
	{
		checkArgs(args);
		if(writeFN==null)
			writeFN=FileUtils.replaceEnd(countsFN, threshold+".txt.gz");
		
		HashMap<String, Double> mappingPercentage = FileUtils.readDoublehash(mappingPercentageFN);
		BufferedReader reader = FileUtils.createReader(mappingPercentageFN);
		String line = null;
		while((line=reader.readLine())!=null)
		{
			String[] cells = line.split("\t");
			String fn = cells[0];
			double mapping = Double.parseDouble(cells[2].replace(",", ""));
			double total = Double.parseDouble(cells[1].replace(",", ""));
			double percentage = mapping/total;
			mappingPercentage.put(fn, percentage);
		}
		MyMatrix counts = new MyMatrix(countsFN);
		HashMap<String, Double> mappingAboveCutoff = getN_AboveCutoff(mappingPercentage,counts);
		
		MyMatrix output = getMappingAboveCutoff(mappingAboveCutoff,counts);
		
		output.write(writeFN);
		System.out.println("mappingPercentageFile file written to: " + writeFN );
	}
	private static MyMatrix getMappingAboveCutoff(HashMap<String, Double> mappingAboveCutoff, MyMatrix counts) {
		
		MyMatrix output = new MyMatrix(counts.rows(),mappingAboveCutoff.size());
		output.rowNames=counts.getRowHeaders();

		int outCol = 0;
		for(int c = 0; c < counts.cols(); c++)
		{
			if(!mappingAboveCutoff.containsKey(counts.colNames[c]))
				continue;
			
			for(int r = 0; r < counts.rows(); r++)
			{
				output.values[r][outCol]=counts.values[r][c];
			}
			output.colNames[outCol]=counts.getColHeaders()[c];
			outCol++;
		}
		return output;
	}
	private static HashMap<String, Double> getN_AboveCutoff(HashMap<String, Double> mappingPercentage, MyMatrix counts) {
		Set<String> fileNames = mappingPercentage.keySet();
		HashMap<String, Double> mappingAboveCutoff = new HashMap<String, Double>();
		for(String errorFN : fileNames)
		{
			File errorFile = new File(errorFN);
			File path = new File(errorFile.getParent());
			String colName =  path.getName();
			double percentage = mappingPercentage.get(errorFN);
			if(percentage>threshold && counts.getColHash().containsKey(colName))
				mappingAboveCutoff.put(colName, percentage);
		}
		return mappingAboveCutoff;
	}
	static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko") && args.length < 1)
			return;
		if(args.length < 1)
		{
			System.out.println("Script requires the following argumetns:\n"
					+ "1. countsFN=<countfilename.txt> - File with the counts for every sample\n"
					+ "2. mappingPercentageFN=<mappingPercentageFN.txt> - File with .err filenames(incl directory) in 1st column and mapping percentage in 2nd column\n"
					+ "3. writeFN=<writeFN.txt> - Filename of the output file\n"
					+ "4. threshold=<0.7> - Value indicating the threshold <default=0.7>\n");
			System.exit(1);
		}
		
		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			System.out.println(args[a]);
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase())
			{
//				var = new JSONutil<Vars>().read(var.JSON_FN, var);
				case "countsFN":
					countsFN =value;
				break;
				case "filename":
					countsFN =value;
					break;
				case "mappingpercentagefn":
					mappingPercentageFN =value;
					break;
				case "writefn":
					writeFN =value;
					break;
				case "threshold":
					threshold = Double.parseDouble(value);
					break;
				default:
					System.out.println("Incorrect argument supplied:\n"+ args[a] +"\nexiting");
					System.exit(1);
			}
		}
	}
}
