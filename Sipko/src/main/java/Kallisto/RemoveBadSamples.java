package Kallisto;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Hashtable;

import MatrixScripts.MyMatrix;
import Tools.FileUtils;

import java.io.FileWriter;

public class RemoveBadSamples {
	//removes samples with expression for less then 10% of the genes
	
	public static String expressionFN = "E:/Groningen/Data/RNAseq_clinic/DownSamples/GRCh38/31.07.pc1.illumina.genes.expressed_DownSamples_ScaffoldGenesRemoved.txt";
	public static String writeFolder = null;
	public static int cutOff =0;//number of reads/transcripts that should be expressed (if 0, will become based on minPercentageExpressed)
	public static int minPercentageExpressed= 10;
	
	public static void main (String[] args) throws IOException
	{
		checkArgs(args);
		if(writeFolder == null)
			writeFolder = expressionFN.replace(".txt", "").replace(".gz", "");
		MyMatrix expression = new MyMatrix(expressionFN);
		
		Hashtable<String,Integer> genesExpressed = new Hashtable<String,Integer>();
		expression.putGenesOnRows();
		System.out.println("rows = " + expression.rows());
		if(cutOff == 0)
			cutOff = (int)(expression.rows()/(100/((double)minPercentageExpressed)));
		System.out.println("cutOff = " + cutOff);
		int nGood = 0;
		
		//determine which samples are good
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(FileUtils.replaceEnd(expressionFN, "_genesExpressed.txt"))));
		writer.write("Sample\tExpressed\tTotal\n");
		for(int c = 0; c < expression.cols(); c++)
		{
			int expressed = 0;
			for(int r = 0; r < expression.rows(); r++)
			{
				if(expression.matrix.get(r,c) > 0)
					expressed++;
			}
			writer.write(expression.getColHeaders()[c]+"\t"+expressed+"\t"+expression.rows()+"\n");
			genesExpressed.put(expression.getColHeaders()[c], expressed);
			if(expressed > cutOff && !expression.colNames[c].contains("DISCARDED"))
				nGood++;
			else
				System.out.println("Sample removed:" + expression.getColHeaders()[c]);//samples that do not make the cutoff (e.g. the bad samples)
		}
		System.out.println("nGood Samples ="+ nGood + " nRemoved = "+ (expression.cols()-nGood));
		writer.close();
		
		writeLeftOverMatrix(expression, nGood, genesExpressed);

		System.out.println("Done, file written to: " + writeFolder);
	}
	
	private static MyMatrix writeLeftOverMatrix(MyMatrix expression, int nGood, Hashtable<String, Integer> genesExpressed) throws IOException {
		MyMatrix result = new MyMatrix(expression.rows(),nGood);
		result.setRowHeaders(expression.getRowHeaders());
		
		int outCol = 0;
		MyMatrix nExpressed = new MyMatrix(expression.rows(),2);
		nExpressed.rowNames=expression.getRowHeaders();
		nExpressed.colNames=new String[]{"expressedN","total"};
		for(int c = 0; c < expression.cols(); c++)
		{
			if(genesExpressed.get(expression.getColHeaders()[c]) <= cutOff || expression.colNames[c].contains("DISCARDED"))
				continue;
			result.setColHeader(outCol,expression.getColHeaders()[c]);
			for(int r = 0; r < result.rows(); r++)
			{
				double value = expression.matrix.get(r, c);
				result.matrix.set(r, outCol, value);
				if(value >0)
					nExpressed.values[r][0]++;
				if(c==0)
					nExpressed.values[r][1]=nGood;
			}
			outCol++;
		}
		nExpressed.write(writeFolder+"_genesExpressedOutput.txt");
		File exp = new File(expressionFN);
		String newName = exp.getName().replace(".txt", "").replace(".gz", "")+"_"+minPercentageExpressed+"PercentTranscriptsExpressed.txt.gz";
		result.write(writeFolder+newName);
		return result;
	}

	static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko") && args.length < 1)
			return;
		if(args.length < 1)
		{
			System.out.println("Script requires the following argumetns:\n"
					+ "1. filename=<expressionFN.txt> - Expression file with genes on rows samples on columns\n"
					+ "2. writeFolder=<writeFolderFN.txt> - Folder where the files will be written (default=filename+_removedBadSamples.txt.gz)\n"
					+ "3. minPercentageExpressed=<10> - Percentage of transcripts/genes that need to be expressed for sample to be included (default=10)\n");
			System.exit(1);
		}
		
		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase())
			{
				case "countFN":
					expressionFN =value;
				break;
				case "filename":
					expressionFN =value;
					break;
				case "writefolder":
					writeFolder = value;
					break;
				case "minpercentageexpressed":
					minPercentageExpressed = Integer.parseInt(value);
					break;
				default:
					System.out.println("Incorrect argument supplied:\n"+ args[a] +"\nexiting");
					System.exit(1);
			}
		}
	}
}
