package STAR;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.Writer;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashMap;
import java.util.function.Consumer;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import Tools.FileUtils;
import Tools.JSONutil;
import Tools.Script;

public class SpliceRatioCalculator extends Script<SpliceRatioCalculator>
{
	//Normalizes the number of reads overlapping a splice variant for the total number of reads overlapping the entire gene
	//file with the read counts needs to be in the same folder as the splice file
	/**
	 * 
	 */
	private static final long serialVersionUID = -6399070446175447150L;
	private String annotationFNComment = "/root/directory/annotation.txt; MANDATORY // A file containing ensembl IDs in the 1st column, chromosome in 2nd, start position in 3rd, end position in 4th column";
	private String annotationFN = null;
	private String ensgToGeneSymbolFNComment = "/root/directory/ensemblToGeneSymbolFN.txt; MANDATORY // A file containing ensembl IDs in the first column and corresponding Gene Symbols in the second column";
	private String ensgToGeneSymbolFN = null;//"E:/Groningen/Data/Annotation/GRCh38/EnsgToGeneSymbol.txt";//"/groups/umcg-gdio/tmp04/umcg-svandam/Data/RNAseq/Annotation/GCRh38/EnsgToGeneSymbol.txt";//
	
	public SpliceRatioCalculator() {
	}

	public String run(String spliceFile, HashMap<String,String[]> spliceSiteToGene) throws FileNotFoundException, IOException
	{
		return run(spliceFile, spliceSiteToGene, null);
	}
	
	public String run(String spliceFile, HashMap<String,String[]> spliceSiteToGene, String geneExpressionFN) throws FileNotFoundException, IOException
	{
		String writeFN = FileUtils.addBeforeExtention(spliceFile, "_GeneSymbols");
		if(geneExpressionFN == null)
			geneExpressionFN = spliceFile.replace(".txt", "_spliceReadsPerGene.txt").replace(".tab", "_spliceReadsPerGene.tab");
		
		addGeneNamesAndCountReadsPerGene(spliceFile,geneExpressionFN, spliceSiteToGene);
		
		BufferedReader spliceReader = FileUtils.createReader(writeFN);
		String correctedWriteFN = writeFN.replace(".tab", "_RatiosPerGene.tab");
		BufferedWriter spliceCorrectedWriter = FileUtils.createWriter(correctedWriteFN);
		//make hash that contains the number of reads per gene
		HashMap<String, Double> geneToReads = FileUtils.makeHash(geneExpressionFN,1);
		
		spliceReader.lines().forEach(line -> calculateRatiosPerGene(line, spliceCorrectedWriter, geneToReads));
		spliceReader.close();
		spliceCorrectedWriter.close();
		System.out.println("File written to: " + correctedWriteFN);
		return correctedWriteFN;
	}

	private void calculateRatiosPerGene(String line, BufferedWriter spliceCorrectedWriter, HashMap<String, Double> geneToReads) 
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
	
	
	public void addGeneNamesAndCountReadsPerGene(String spliceFN, 
			String readsPerGeneWriteFN, 
			HashMap<String, String[]> spliceSiteToGene)
	{
		addGeneNamesAndCountReadsPerGene(spliceFN, readsPerGeneWriteFN, spliceSiteToGene, null, true); 
	}
	
	public void addGeneNamesAndCountReadsPerGene(String spliceFN, 
			String readsPerGeneWriteFN, 
			HashMap<String, String[]> spliceSiteToGene, 
			String writeFN, 
			boolean afterLine) 
	{
		try 
		{ 
			HashMap<String,Integer> readsPerGene = new HashMap<String,Integer>();
			HashMap<String,ArrayList<Gene>> genome = new HashMap<String,ArrayList<Gene>>(); 
			BufferedReader annotationReader = FileUtils.createReader(getAnnotationFN()); 
			HashMap<String,String> ensemblToGeneSymbol = FileUtils.readStringStringHash(getEnsgToGeneSymbolFN());
			//add all genes to the genome
			//System.out.println("Building genome reference set");
			annotationReader.lines().skip(1).forEach(line -> addToGenome(line,genome,ensemblToGeneSymbol));
			//add genes to the spliceFile
			System.out.println("Adding genes to the spliceFile");
			if(writeFN == null)
				writeFN = FileUtils.addBeforeExtention(spliceFN,"_GeneSymbols");
			if(readsPerGeneWriteFN == null)
				readsPerGeneWriteFN = spliceFN.replace(".txt", "_SpliceReadsPerGene.txt").replace(".tab", "_SpliceReadsPerGene.tab");
			
			BufferedReader spliceReader = FileUtils.createReader(spliceFN);
			BufferedWriter spliceWriter = FileUtils.createWriter(writeFN);
			//add the header
			String firstline = spliceReader.readLine();
			if(firstline.contains("First base of the intron (1-based)"))
			{
				if(!afterLine)
					spliceWriter.write("EnsemblIDs\tGeneSymbols\t"+firstline+"\n");
				else
					spliceWriter.write(firstline+"\t"+"EnsemblIDs\tGeneSymbols"+"\n");
			}
			else
				spliceReader = FileUtils.createReader(spliceFN);//start the stream at the first line again
			
			spliceReader.lines().forEach(line -> 
			{
				String[] cells = line.split("\t");
				String chromLoc = cells[0]+"\t"+cells[1]+"\t"+cells[2]+"\t"+cells[3]+"\t"+cells[4];
				String chromosomeName = cells[0];
				addGeneNames_And_CountReadsPerGene(line, cells, chromosomeName, chromLoc, genome,spliceWriter,readsPerGene, spliceSiteToGene, afterLine);
			});

			BufferedWriter readsPerGeneWriter = FileUtils.createWriter(readsPerGeneWriteFN);
			readsPerGene.forEach((gene,reads)-> writeReadsPerGene(gene,reads,readsPerGeneWriter));
			System.out.println("ReadsPerGeneWrite written to: " + readsPerGeneWriteFN);
			readsPerGeneWriter.close();
			spliceReader.close();
			spliceWriter.close();
			
		} catch (Exception e) {	e.printStackTrace();}	
	}
	
	public void writeReadsPerGene(String gene, int reads, BufferedWriter readsPerGeneWriter)
	{
		try {
			readsPerGeneWriter.write(gene+"\t"+reads+"\n");
		} catch (IOException e) {e.printStackTrace();}
	}
	
	
	private void addGeneNames_And_CountReadsPerGene(String line, 
													String[] cells, 
													String chromosomeName, 
													String chromLoc, 
													HashMap<String, ArrayList<Gene>> genome,
													BufferedWriter spliceWriter, 
													HashMap<String, Integer> readsPerGene, 
													HashMap<String,String[]> spliceSiteToGene, 
													boolean afterLine) 
	{
		try 
		{
			int reads = Integer.parseInt(cells[6]);
			int start,end;
			
			if(genome.containsKey(chromosomeName))
			{
				ArrayList<Gene> chromosome = genome.get(chromosomeName);
				String[] genesIDsToAdd = spliceSiteToGene.get(chromLoc);
				
				if(genesIDsToAdd == null)
				{
					start = Integer.parseInt(cells[1]);
					end = Integer.parseInt(cells[2]);
					String[] getGenes= new String[2];
					chromosome.forEach(gene -> addGeneString(gene,start,end, getGenes));//puts results into getGenes
					genesIDsToAdd=getGenes;
					spliceSiteToGene.put(chromLoc, genesIDsToAdd);
				}

				if(!afterLine)
					spliceWriter.write(genesIDsToAdd[0]+"\t"+genesIDsToAdd[1]+"\t"+line+"\n");
				else
					spliceWriter.write(line+"\t"+genesIDsToAdd[0]+"\t"+genesIDsToAdd[1]+"\n");
				
				//count the number of spliced reads per gene (used to calculate the ratios of splice variants in a different script)
				if(genesIDsToAdd[0]!=null)
				{
					String[] geneIDs= genesIDsToAdd[0].split(",");
					Stream.of(geneIDs).forEach(geneID -> 
					{
						if(readsPerGene.containsKey(geneID))
						{
							int number = readsPerGene.get(geneID);
							readsPerGene.put(geneID, number+reads);
						}
						else
							readsPerGene.put(geneID, reads);
					});
				}
			}
			else
			{
				if(!afterLine)
					spliceWriter.write("\t-\t-"+line+"\n");
				else
					spliceWriter.write(line+"\t-\t-"+"\n");
				spliceSiteToGene.put(chromLoc, new String[]{"-","-"});
			}
		} catch (Exception e) {e.printStackTrace();}
	}


	private String[] addGeneString(Gene gene, int start, int end, String[] genesIDsToAdd) 
	{
		if(gene.inGene(start) || gene.inGene(end))
		{
			if(genesIDsToAdd[0]==null)
			{
				genesIDsToAdd[0]=gene.ensemblID;
				genesIDsToAdd[1]=gene.geneSymbol;
			}
			else
			{
				genesIDsToAdd[0]+=","+gene.ensemblID;
				genesIDsToAdd[1]+=","+gene.geneSymbol;
			}
		}
		return genesIDsToAdd;
	}

	private void addToGenome(String line, HashMap<String, ArrayList<Gene>> genome, HashMap<String,String> ensemblToGeneSymbol)
	{
		String[] cells = line.split("\t");
		String ensemblID = cells[0];
		String geneSymbol = ensemblToGeneSymbol.get(ensemblID);
		String chromosomeName= cells[1];
		int start = Integer.parseInt(cells[2]);
		int end = Integer.parseInt(cells[3]);
		
		Gene gene = new Gene(ensemblID, geneSymbol, chromosomeName, start, end);
		if(genome.containsKey(chromosomeName))
		{
			ArrayList<Gene> chromosome = genome.get(chromosomeName);
			chromosome.add(gene);
			genome.put(chromosomeName, chromosome);
		}
		else
		{
			ArrayList<Gene> chromosome = new ArrayList<Gene>();
			chromosome.add(gene);
			genome.put(chromosomeName, chromosome);
		}
		
	}
	
	private void writeNew(String cell, BufferedWriter spliceCorrectedWriter) {
		try {
			spliceCorrectedWriter.write(cell+"\t");
		} catch (IOException e) {e.printStackTrace();}
	}

	public String getAnnotationFN() {
		return annotationFN;
	}

	public void setAnnotationFN(String annotationFN) {
		this.annotationFN = annotationFN;
	}

	public String getEnsgToGeneSymbolFN() {
		return ensgToGeneSymbolFN;
	}

	public void setEnsgToGeneSymbolFN(String ensgToGeneSymbolFN) {
		this.ensgToGeneSymbolFN = ensgToGeneSymbolFN;
	}
}
