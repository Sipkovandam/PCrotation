package Analyses;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

import Tools.FileUtils;

public class AveragesPerRow {
	// calculates the averages and standard deviations per row, excluding NaN and INFINITE numbers.
	
	static String fileName = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/TESTexpression.txt";
	//static String writeFN = "/Volumes/Promise_RAID/GeneNetwork/Sipko/04-2016/22000Samples_rLog_1.0_Correl/gene_correlation_absoluteAverages.txt";
	static String writeFN = null;
	static boolean absolute = false;

	public static void main(String[] args) throws IOException
	{
		
//		String fileName = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/TESTexpression/SAMPLE_covariance.txt";
//		String writeFN = fileName.replace(".txt", "_absoluteAvg.txt");
		checkArgs(args);
		if(writeFN == null)
			writeFN = fileName.replace(".txt", "").replace(".gz", "")+"_rowAverages.txt";
		
		BufferedReader reader = FileUtils.createReader(fileName);
		BufferedWriter writerAverages = new BufferedWriter(new FileWriter(new File(writeFN)));
		BufferedWriter writerSTdevs = new BufferedWriter(new FileWriter(new File(writeFN.replace("_rowAverages.txt", "row_STDevs.txt"))));
		
		String line = reader.readLine();//skip header col
		writerAverages.write("RowName\tAverage\n");
		writerSTdevs.write("RowName\tStandard deviation\n");
		
		while((line = reader.readLine()) !=null)
		{
			String[] eles = line.split("\t");
			double total = 0;
			double totalEles = (eles.length-1);

			int nValid = 0;
			for(int e = 1; e < eles.length; e++)// this makes this script a lot slower having to parse the same number 2X, but w/e
			{
				double number = Double.parseDouble(eles[e]);
				if(!Double.isNaN(number) && !Double.isInfinite(number))
					nValid++;
			}
			
			double[] rowValues = new double[nValid];//used to calculate stdevs
			int index = 0;
			
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
				
				rowValues[index] = number;
				index++;
			}
			double average = total/totalEles;
			double variance = org.apache.commons.math3.stat.StatUtils.variance(rowValues);
			double standardDev = java.lang.Math.pow(variance,0.5);
			writerAverages.write(eles[0]+"\t" + average+"\n");
			writerSTdevs.write(eles[0]+"\t" + standardDev +"\n");
		}
		writerAverages.close();
		writerSTdevs.close();
		reader.close();
		
		
		System.out.println("Done! File written to: " + writeFN);
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
					System.out.println("Incorrect argument supplied: " + args[a] + " exiting");
					System.exit(1);
			}
		}
		if(args.length < 1)
		{
			System.out.println("Wrong arguments, requires:\n"
					+ "1. fileName=<filename.txt> - name of the file.\n"
					+ "2. writeFN=<writeFileName.txt> - name of the output file (default=<fileName>_rowAverages.txt\n"
					+ "3. absolute=<true/false> - if true values will be made absolute before taking average");
			System.exit(1);
		}

	}
}
