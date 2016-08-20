package PrepData;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import PCA.MatrixStruct;

public class GetSamplesWithEmptyCells 
{
	static String fn = "E:/Groningen/Data/Juha/Calculon/JuhaMerged/test/TestAwk.txt";
	static String writeFN = "E:/Groningen/Data/Juha/Calculon/JuhaMerged/test/TestJava.txt";
	
	public static void main (String[] args) throws IOException
	{	
		checkArgs(args);
		File file = new File(fn);
		BufferedReader reader = new BufferedReader(new FileReader(file));
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(writeFN)));
		String line = reader.readLine();//header line
		
		int length = line.split("\t").length;
		ArrayList<String> failedSamples = new ArrayList<String>();
		failedSamples.add(line);
		while((line=reader.readLine())!= null)
		{
			String[] values = line.split("\t");
			if(values.length != length)
			{
				writer.write(values[0]+"\n");
				failedSamples.add(line);
				continue;
			}
				
			for(int v = 1; v < values.length; v++)
				if(values[v].length()==0)
				{
					writer.write(values[0]+"\n");
					failedSamples.add(line);
					break;
				}
		}
		reader.close();
		writer.close();
		
		int longest = 0;
		for(int l = 0; l < failedSamples.size();l++)
		{
			String sample = failedSamples.get(l);
			String[] cells = sample.split("\t",-1);
			if(cells.length>longest)
				longest = cells.length+1;
		}
		//make an array where the columns are the samples (easier to read in excel)
		String[][] samples = new String[longest][failedSamples.size()];
		for(int l = 0; l < failedSamples.size();l++)
		{
			String sample = failedSamples.get(l);
			String[] cells = sample.split("\t",-1);
			for(int c = 0; c< cells.length; c++)
			{
				samples[c][l] = cells[c];
				if(cells.length>samples.length)
					System.out.println("cellsLength =" + cells.length + " samples.length=" + samples.length);
			}
		}
		
		BufferedWriter sampleWriter = new BufferedWriter(new FileWriter(new File (writeFN.replace(".txt", "_sampleValues.txt"))));
		for(int r = 0; r < samples.length; r++)
		{
			if(samples[r][0] != null)
				sampleWriter.write(samples[r][0]);
			for(int c = 1; c < samples[0].length; c++)
				if(samples[r][0] != null)
					sampleWriter.write("\t"+samples[r][c]);
				else
					sampleWriter.write("\t");
			sampleWriter.write("\n");
		}
		
		
		sampleWriter.close();
		System.out.println("File written to: " + writeFN);
	}
	static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
			return;
		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase()){
				case "filename":
					fn =value;
					break;
				case "writefn":
					writeFN = value;
					break;
				default:
					checkArgs(args);
					System.out.println("Incorrect argument supplied: " + args[a] + " exiting");
					System.exit(1);
			}
		}
		if(args.length < 1)
		{
			System.out.println("Wrong arguments, requires:\n"
					+ "1. fileName=<filename.txt> - name of the file.\n"
					+ "2. writeFN=<writeFileName.txt> - name of the output file (default=<fileName>_rowAverages.txt\n");
			System.exit(1);
		}

	}
}
