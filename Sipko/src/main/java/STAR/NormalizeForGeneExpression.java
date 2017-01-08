package STAR;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.Writer;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.function.Consumer;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import Tools.FileUtils;
import Tools.JSONutil;

public class NormalizeForGeneExpression 
{
	//Normalizes the number of reads overlapping a splice variant for the total number of reads overlapping the entire gene
	//file with the read counts needs to be in the same folder as the splice file
	
	public static void main(String[] args) throws FileNotFoundException, IOException 
	{
		String spliceFile =  "E:/Groningen/Test/SpliceNormalization/SJ.out.tab";
		String expressionFN = null;//"E:/Groningen/Test/SpliceNormalization/SJ_readsPerGene.tab.out.tab";
		
		Variables v = checkArgs(args);
//		if(var.writeFolder ==null)
//			var.writeFolder = new File(var.writeFN).getParent()+"/";

		HashMap<String, String[]> startEndToGene = new HashMap<String,String[]>();
		run(v, spliceFile, startEndToGene, expressionFN);
	}
	
	public static String run(Variables v, String spliceFile, HashMap<String,String[]> spliceSiteToGene) throws FileNotFoundException, IOException
	{
		return run(v, spliceFile, spliceSiteToGene, null);
	}
	
	public static String run(Variables v, String spliceFile, HashMap<String,String[]> spliceSiteToGene, String geneExpressionFN) throws FileNotFoundException, IOException
	{
		String writeFN = spliceFile.replace(".txt", "_GeneSymbols.txt").replace(".tab", "_GeneSymbols.tab");
		if(geneExpressionFN == null)
			geneExpressionFN = spliceFile.replace(".txt", "_spliceReadsPerGene.txt").replace(".tab", "_spliceReadsPerGene.tab");
		
		SpliceSitesPerGene.addGeneNamesAndCountReadsPerSplice(v, spliceFile,geneExpressionFN, spliceSiteToGene);
		
//		SpliceSitesPerGene.main(new String[]{ "splicefn="+spliceFile
//				, "annotationfn="+v.annotationFN 
//				, "afterline=true"
//				, "writeFN="+writeFN
//				, "readsPerGeneWriteFN="+expressionFN
//				}, startEndToGene);
		
		BufferedReader spliceReader = FileUtils.createReader(writeFN);
		String correctedWriteFN = writeFN.replace(".tab", "_Normalized.tab");
		BufferedWriter spliceCorrectedWriter = FileUtils.createWriter(correctedWriteFN);
		//make hash that contains the number of reads per gene
		HashMap<String, Double> geneToReads = FileUtils.makeHash(geneExpressionFN,1);
		
		spliceReader.lines().forEach(line -> normalize(line, spliceCorrectedWriter, geneToReads));
		spliceReader.close();
		spliceCorrectedWriter.close();
		System.out.println("File written to: " + correctedWriteFN);
		return correctedWriteFN;
	}

	private static void normalize(String line, BufferedWriter spliceCorrectedWriter, HashMap<String, Double> geneToReads) 
	{
		String[] cells = line.split("\t");
		Double spliceCount = Double.parseDouble(cells[6]);
		String gene = "DoesNotExist";
		if(cells.length>8)
			gene=cells[9];
		
		Stream.of(cells).forEach(cell->writeNew(cell,spliceCorrectedWriter));
		
		String[] genes = gene.split(",");
		try 
		{
			double geneCounts = 0;
			for(int g = 0; g < genes.length; g++)
			{
				if(geneToReads.containsKey(genes[g]))
					geneCounts+= geneToReads.get(genes[g]);
	
				if(g==0)
					spliceCorrectedWriter.write(""+spliceCount/geneCounts);
				else
					spliceCorrectedWriter.write(","+spliceCount/geneCounts);
	//				DecimalFormat format = new DecimalFormat("#.####");
	//				spliceCorrectedWriter.write(format.format(spliceCount/geneCounts)+"\n");
				
			}
			spliceCorrectedWriter.write("\n");
		} catch (IOException e) {e.printStackTrace();};
	}

	private static void writeNew(String cell, BufferedWriter spliceCorrectedWriter) {
		try {
			spliceCorrectedWriter.write(cell+"\t");
		} catch (IOException e) {e.printStackTrace();}
	}

	static Variables checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko") && args.length < 1)
			return new Variables();
		if(args.length < 1)
		{
			System.out.println("Script requires the following argumetns:\n"
					+ "1. filename=<expressionFN.txt> - Expression file with genes on rows samples on columns\n");
			System.exit(1);
		}
		
		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase())
			{
				case "json":
					return Variables.readVars(value);
				default:
					System.out.println("Incorrect argument supplied:\n"+ args[a] +"\nexiting");
					System.exit(1);
			}
		}
		return null;
	}
}
