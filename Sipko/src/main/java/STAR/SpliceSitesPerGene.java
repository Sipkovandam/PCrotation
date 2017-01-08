package STAR;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.stream.Stream;

import JuhaPCA.FileUtil;
import PCA.Matrix;
import Tools.FileUtils;
import Tools.JSONutil;

public class SpliceSitesPerGene 
{
	//Adds gene names in two columns to the splice file
	
	public static void main(String[] args) throws FileNotFoundException, IOException
	{
		String spliceFile = "E:/Groningen/Splicing/Results/L1CAM_CACGAT/SJ_zScoresAdded.out.tab";
		String readsPerGeneWriteFN = null;
		String writeFN = spliceFile.replace(".txt", "_GeneSymbols.txt").replace(".tab", "_GeneSymbols.tab");
		boolean afterLine = true;
		checkArgs(args);
		HashMap<String,String[]> startEndToGene= new HashMap<String,String[]>();
		Variables v = new Variables();
		addGeneNamesAndCountReadsPerSplice(v, spliceFile, readsPerGeneWriteFN, startEndToGene, writeFN, afterLine);
		
		System.out.println("Splice file, with gene names added, written to: " + writeFN);
	}
	
	public static void addGeneNamesAndCountReadsPerSplice(Variables v, 
			  String spliceFN, 
			  String readsPerGeneWriteFN, 
			  HashMap<String, String[]> spliceSiteToGene)
	{
		addGeneNamesAndCountReadsPerSplice(v, spliceFN, readsPerGeneWriteFN, spliceSiteToGene, null, true); 
	}
	
	public static void addGeneNamesAndCountReadsPerSplice(Variables v, 
														  String spliceFN, 
														  String readsPerGeneWriteFN, 
														  HashMap<String, String[]> spliceSiteToGene, 
														  String writeFN, 
														  boolean afterLine) 
	{
		try 
		{ 
			HashMap<String,Integer> readsPerGene = new HashMap<String,Integer>();
			Hashtable<String,ArrayList<Gene>> genome = new Hashtable<String,ArrayList<Gene>>(); 
			BufferedReader annotationReader = FileUtils.createReader(v.getAnnotationFN()); 
			Hashtable<String,String> ensemblToGeneSymbol = FileUtils.readStringStringHash(v.getEnsgToGeneSymbolFN());
			//add all genes to the genome
			//System.out.println("Building genome reference set");
			annotationReader.lines().skip(1).forEach(line -> addToGenome(line,genome,ensemblToGeneSymbol));
			//add genes to the spliceFile
			System.out.println("Adding genes to the spliceFile");
			if(writeFN == null)
				writeFN = spliceFN.replace(".txt", "_GeneSymbols.txt").replace(".tab", "_GeneSymbols.tab");
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
				addGeneNames_And_CountReadsPerGene(v, line, cells, chromosomeName, chromLoc, genome,spliceWriter,readsPerGene, spliceSiteToGene, afterLine);
			});

			BufferedWriter readsPerGeneWriter = FileUtils.createWriter(readsPerGeneWriteFN);
			readsPerGene.forEach((gene,reads)-> writeReadsPerGene(gene,reads,readsPerGeneWriter));
			System.out.println("ReadsPerGeneWrite written to: " + readsPerGeneWriteFN);
			readsPerGeneWriter.close();
			spliceReader.close();
			spliceWriter.close();
			
		} catch (Exception e) {	e.printStackTrace();}	
	}

	public static void writeReadsPerGene(String gene, int reads, BufferedWriter readsPerGeneWriter)
	{
		try {
			readsPerGeneWriter.write(gene+"\t"+reads+"\n");
		} catch (IOException e) {e.printStackTrace();}
	}
	
	private static void addGeneNames_And_CountReadsPerGene(Variables v, 
														  String line, 
														  String[] cells, 
														  String chromosomeName, 
														  String chromLoc, 
														  Hashtable<String, ArrayList<Gene>> genome,
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
							//return readsPerGene.get(genesIDsToAdd[0])+reads;
							int number = readsPerGene.get(geneID);
							readsPerGene.put(geneID, number+reads);
						}
						else
							//return reads;
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

	private static String[] addGeneString(Gene gene, int start, int end, String[] genesIDsToAdd) 
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

	private static void addToGenome(String line, Hashtable<String, ArrayList<Gene>> genome, Hashtable<String,String> ensemblToGeneSymbol)
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

	static Variables checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko") && args.length < 1)
			return new Variables();
		if(args.length < 1)
		{
			System.out.println("Script requires the following argumetns:\n"
					+ "1. filename=<expressionFN.txt> - Expression file with genes on rows samples on columns\n"
					+ "2. writeFolder=<writeFolderFN.txt> - Folder where the files will be written (default=parentFolder(input.txt))\n"
					+ "3. geoFN=<geoFn.txt> - Optional file with geometric mean per gene to use (default=null)\n");
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
//				case "splicefn":
//					var.spliceFN = value;
//					break;
//				case "splicefile":
//					var.spliceFN = value;
//					break;
//				case "writefn":
//					var.writeFN = value;
//					break;
//				case "annotationfn":
//					var.annotationFN = value;
//					break;
//				case "afterline":
//					var.afterLine = Boolean.parseBoolean(value);
//					break;
//				case "readspergenewritefn":
//					var.readsPerGeneWriteFN = value;
//					System.out.println("test=" + var.readsPerGeneWriteFN );
//					break;
				default:
					System.out.println("Incorrect argument supplied:\n"+ args[a] +"\nexiting");
					System.exit(1);
			}
		}
		return null;
	}
}
