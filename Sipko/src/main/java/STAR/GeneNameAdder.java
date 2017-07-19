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

public class GeneNameAdder extends Script<GeneNameAdder>
{
	//Normalizes the number of reads overlapping a splice variant for the total number of reads overlapping the entire gene
	//file with the read counts needs to be in the same folder as the splice file
	/**
	 * 
	 */
	private static final long serialVersionUID = -6399070446175447150L;
	private String spliceFnComment = "/root/directory/SpliceJunction_expression.txt; MANDATORY // A file containing the expression per splice junction. Works with 2 formats";
	private String spliceFn = null;//"E:/Groningen/Data/Annotation/GRCh38/EnsgToGeneSymbol.txt";//"/groups/umcg-gdio/tmp04/umcg-svandam/Data/RNAseq/Annotation/GCRh38/EnsgToGeneSymbol.txt";//
	private String annotationFnComment = "/root/directory/annotation.txt; MANDATORY // A file containing ensembl IDs in the 1st column, chromosome in 2nd, start position in 3rd, end position in 4th column";
	private String annotationFn = null;
	private String ensgToGeneSymbolFNComment = "/root/directory/ensemblToGeneSymbolFN.txt; MANDATORY // A file containing ensembl IDs in the first column and corresponding Gene Symbols in the second column";
	private String ensgToGeneSymbolFN = null;//"E:/Groningen/Data/Annotation/GRCh38/EnsgToGeneSymbol.txt";//"/groups/umcg-gdio/tmp04/umcg-svandam/Data/RNAseq/Annotation/GCRh38/EnsgToGeneSymbol.txt";//
	private String writeFn = null;
	
	public GeneNameAdder()
	{
	}
	public void run()
	{
		if(writeFn == null)
			writeFn = FileUtils.addBeforeExtention(	spliceFn,
				"_geneNamesAdded");
		addGeneNames(	spliceFn,
		new HashMap<String, String[]>(),
		writeFn,
		false);
	}              

	public String run(	String spliceFile,
						HashMap<String, String[]> spliceSiteToGene) throws FileNotFoundException, IOException
	{
		String writeFN = FileUtils.addBeforeExtention(	spliceFile,
														"_GeneSymbols");
		addGeneNames(	spliceFile,
						spliceSiteToGene,
						writeFN);
		return writeFN;
	}

	public void addGeneNames(	String spliceFN,
								HashMap<String, String[]> spliceSiteToGene,
								String writeFN)
	{
		addGeneNames(	spliceFN,
						spliceSiteToGene,
						writeFN,
						true);
	}

	public void addGeneNames(	String spliceFN,
								HashMap<String, String[]> spliceSiteToGene,
								String writeFN,
								boolean afterLine)
	{
		try
		{
			//Holds all chromosomes. For each chromosome has an arraylist of all Genes present on the chromosome.
			HashMap<String, ArrayList<Gene>> genome = new HashMap<String, ArrayList<Gene>>();
			BufferedReader annotationReader = FileUtils.createReader(getAnnotationFN());
			HashMap<String, String> ensemblToGeneSymbol = FileUtils.readStringStringHash(getEnsgToGeneSymbolFN());
			
			//add all genes to the genome
			annotationReader.lines().skip(1).forEach(line -> addToGenome(	line,
																			genome,
																			ensemblToGeneSymbol));
			//add genes to the spliceFile
			//System.out.println("Adding genes to the spliceFile");
			BufferedReader spliceReader = FileUtils.createReader(spliceFN);
			BufferedWriter spliceWriter = FileUtils.createWriter(writeFN);
			//add the header
			String firstline = spliceReader.readLine();
			if(firstline ==null)
			{
				p("File is empty:\t " + spliceFN);
				return;
			}
			
			if (firstline.contains("First base of the intron (1-based)"))
			{
				if (!afterLine)
					spliceWriter.write("EnsemblIDs\tGeneSymbols\t" + firstline + "\n");
				else
					spliceWriter.write(firstline + "\t" + "EnsemblIDs\tGeneSymbols" + "\n");
			}
			else if(!firstline.startsWith("\t"))//start the stream at the first line again (in case this is a file that does not have a header)
				spliceReader = FileUtils.createReader(spliceFN);
			else
				if(!afterLine)
					spliceWriter.write("EnsemblIds\t" + "GeneSymbols\t" + firstline+"\n");
				else
					spliceWriter.write(firstline+"EnsemblIds\t" + "GeneSymbols\t" +"\n");

			spliceReader.lines().forEach(line ->
			{
				String[] cells = line.split("\t|_");//"_" is for files in which the file junctions are described in 1 column with "_" as separator between aspects of the junction
				String chromLoc = cells[0];// + "\t" + cells[1] + "\t" + cells[2] + "\t" + cells[3] + "\t" + cells[4];
				for(int c = 1; c < 6; c++)
				{
					chromLoc=chromLoc.concat("_");
					chromLoc=chromLoc.concat(cells[c]);
				}
				//System.out.println(chromLoc);
				String chromosomeName = cells[0];
				
				//adds the genenames to the splice junction files (writes it to a new file using "splicewriter")
				addGeneNames(	line,
								cells,
								chromosomeName,
								chromLoc,
								genome,
								spliceWriter,
								spliceSiteToGene,
								afterLine);
			});
			spliceReader.close();
			spliceWriter.close();

		} catch (Exception e)
		{
			e.printStackTrace();
		}
		p("GeneNames Merged, file written at:\n"+writeFn);
	};

	/** adds the genenames to the splice junction files (writes it to a new file using "splicewriter")
	 * 
	 * @param line This is the line that was readily present in the splice juntion file to which the gene name is to be added
	 * @param cells This is the line split by "/t"
	 * @param chromosomeName The chromosome on which this gene is located according to the splice junction name
	 * @param chromLoc The location on the chromosome on which this splice junction is located
	 * @param genome The parsed genome effectively holding the information of which genes are located where on which chromosome
	 * @param spliceWriter Writes the new file that also has the gene names
	 * @param spliceSiteToGene A hash that stores all splice sites for which the genes have readily been retrieved by comparing it to the genome
	 * @param afterLine A boolean defining whether the gene name should be added at the start or at the end of the line
	 */
	private void addGeneNames(	String line,
								String[] cells,
								String chromosomeName,
								String chromLoc,
								HashMap<String, ArrayList<Gene>> genome,
								BufferedWriter spliceWriter,
								HashMap<String, String[]> spliceSiteToGene,
								boolean afterLine)
	{
		try
		{
			int start, end;

			if (genome.containsKey(chromosomeName))
			{
				ArrayList<Gene> chromosome = genome.get(chromosomeName);
				String[] genesIDsToAdd = spliceSiteToGene.get(chromLoc);

				//if the genes for this location were not readily retrieved from the genome, figure out which gene(s) this splicejunction overlaps
				if (genesIDsToAdd == null)
				{
					start = Integer.parseInt(cells[1]);
					end = Integer.parseInt(cells[2]);
					
					//getGenes[0] has all the ensembl name the splice junction overlaps, getGenes[1] has all the gene_symbols that the splice juction overlaps in a comma separated String
					String[] getGenes = new String[2];
					
					//puts results into getGenes
					chromosome.forEach(gene -> addGeneString(	gene,
																start,
																end,
																getGenes));
					genesIDsToAdd = getGenes;
					spliceSiteToGene.put(	chromLoc,
											genesIDsToAdd);
				}
				
				if (!afterLine)
					spliceWriter.write(genesIDsToAdd[0] + "\t" + genesIDsToAdd[1] + "\t" + line + "\n");
				else
					spliceWriter.write(line + "\t" + genesIDsToAdd[0] + "\t" + genesIDsToAdd[1] + "\n");
			}
			else
			{
				if (!afterLine)
					spliceWriter.write("\t-\t-" + line + "\n");
				else
					spliceWriter.write(line + "\t-\t-" + "\n");
				spliceSiteToGene.put(	chromLoc,
										new String[] { "-", "-" });
			}
		} catch (Exception e)
		{
			e.printStackTrace();
		}
	}

	/** Adds a geneName if it overlaps the splice junction 
	 * 
	 * @param gene The gene to be compared to the splice juntion location
	 * @param start The genomic start location of the splice junction
	 * @param end The genomic end location of the splice junction
	 * @param genesIDsToAdd This function is called itteratively and these are the genes that were already added in a previous iteration
	 * @return An array containing a comma separated list of all ensembl IDs in the first element and the gene symbols in the second element of the array.
	 */
	
	private String[] addGeneString(	Gene gene,
									int start,
									int end,
									String[] genesIDsToAdd)
	{
		if (gene.inGene(start) || gene.inGene(end))
		{
			if (genesIDsToAdd[0] == null)
			{
				genesIDsToAdd[0] = gene.getEnsemblID();
				genesIDsToAdd[1] = gene.getGeneSymbol(true);
			}
			else
			{
				genesIDsToAdd[0] += "," + gene.getEnsemblID();
				genesIDsToAdd[1] += "," + gene.getGeneSymbol(true);
			}
		}
		return genesIDsToAdd;
	}

	private void addToGenome(	String line,
								HashMap<String, ArrayList<Gene>> genome,
								HashMap<String, String> ensemblToGeneSymbol)
	{
		String[] cells = line.split("\t");
		String ensemblID = cells[0];
		String geneSymbol = ensemblToGeneSymbol.get(ensemblID);
		String chromosomeName = cells[1];
		int start = Integer.parseInt(cells[2]);
		int end = Integer.parseInt(cells[3]);

		Gene gene = new Gene(	ensemblID,
								geneSymbol,
								chromosomeName,
								start,
								end);
		if (genome.containsKey(chromosomeName))
		{
			ArrayList<Gene> chromosome = genome.get(chromosomeName);
			chromosome.add(gene);
			genome.put(	chromosomeName,
						chromosome);
		}
		else
		{
			ArrayList<Gene> chromosome = new ArrayList<Gene>();
			chromosome.add(gene);
			genome.put(	chromosomeName,
						chromosome);
		}

	}

	private void writeNew(	String cell,
							BufferedWriter spliceCorrectedWriter)
	{
		try
		{
			spliceCorrectedWriter.write(cell + "\t");
		} catch (IOException e)
		{
			e.printStackTrace();
		}
	}

	public String getAnnotationFN()
	{
		return annotationFn;
	}

	public void setAnnotationFN(String annotationFN)
	{
		this.annotationFn = annotationFN;
	}

	public String getEnsgToGeneSymbolFN()
	{
		return ensgToGeneSymbolFN;
	}

	public void setEnsgToGeneSymbolFN(String ensgToGeneSymbolFN)
	{
		this.ensgToGeneSymbolFN = ensgToGeneSymbolFN;
	}
	public String getSpliceFn()
	{
		return spliceFn;
	}
	public void setSpliceFn(String spliceFn)
	{
		this.spliceFn = spliceFn;
	}
	public String pcaDir()
	{
		return writeFn;
	}
	public void setWriteFn(String writeFn)
	{
		this.writeFn = writeFn;
	}
	public String getAnnotationFn()
	{
		return annotationFn;
	}
	public void setAnnotationFn(String annotationFn)
	{
		this.annotationFn = annotationFn;
	}
	public String getWriteFn()
	{
		return writeFn;
	}
}
