package Analyses;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class AveragesPerRow {
	
	public static void main(String[] args) throws IOException
	{
		String fileName = "/Volumes/Promise_RAID/GeneNetwork/Sipko/04-2016/22000Samples_rLog_1.0_Correl/gene_correlation.txt";
		String writeFN = "/Volumes/Promise_RAID/GeneNetwork/Sipko/04-2016/22000Samples_rLog_1.0_Correl/gene_correlation_absoluteAverages.txt";
		
//		String fileName = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/TESTexpression/SAMPLE_covariance.txt";
//		String writeFN = fileName.replace(".txt", "_absoluteAvg.txt");
		
		boolean absolute = true;
		
		if(args.length < 1)
			checkArgs(args);
		
		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase()){
				case "filename":
					fileName =value;
					break;
				case "writefn":
					writeFN = value;
					break;
				case "absolute":
					absolute = Boolean.parseBoolean(value);
					break;
				default:
					checkArgs(args);
					System.out.println("Incorrect argument supplied; exiting");
					System.exit(1);
			}
		}
		
		BufferedReader reader = new BufferedReader(new FileReader(new File(fileName)));
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(writeFN)));
		String line = reader.readLine();//skip header col
		writer.write("RowName\tAverage\n");
		while((line = reader.readLine()) !=null)
		{
			String[] eles = line.split("\t");
			double total = 0;
			double totalEles = (eles.length-1);
			for(int e = 1; e < eles.length; e++)
			{
				double number = Double.parseDouble(eles[e]);

				if(Double.isNaN(number) || Double.isInfinite(number))
				{
					totalEles--;
					continue;
				}
				if(absolute)
					total+=Math.abs(number);
				else
					total+=number;
			}
			double average = total/totalEles;
			writer.write(eles[0]+"\t" + average+"\n");
		}
		writer.close();
		reader.close();
		System.out.println("Done! File written to: " + writeFN);
	}
	static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
			return;
		System.out.println("Wrong arguments");
		System.exit(1);
	}
}
