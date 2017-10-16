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
import java.util.HashSet;
import java.util.List;
import java.util.HashMap;
import java.util.function.Consumer;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import Tools.FileUtils;
import Tools.JSONutil;
import Tools.Script;
import umcg.genetica.collections.intervaltree.PerChrIntervalTree;
import umcg.genetica.io.gtf.GffElement;
import umcg.genetica.io.gtf.GtfReader;

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
//	private String annotationFnComment = "/root/directory/annotation.txt; MANDATORY // A file containing ensembl IDs in the 1st column, chromosome in 2nd, start position in 3rd, end position in 4th column";
//	private String annotationFn = null;
//	private String ensgToGeneSymbolFNComment = "/root/directory/ensemblToGeneSymbolFN.txt; MANDATORY // A file containing ensembl IDs in the first column and corresponding Gene Symbols in the second column";
//	private String ensgToGeneSymbolFN = null;//"E:/Groningen/Data/Annotation/GRCh38/EnsgToGeneSymbol.txt";//"/groups/umcg-gdio/tmp04/umcg-svandam/Data/RNAseq/Annotation/GCRh38/EnsgToGeneSymbol.txt";//
	private String writeFn = null;
	private String writeFnNoGeneSymbol = null;
	
	private String gtfFn = null;
	//	private String gtfFn = null;

	public GeneNameAdder()
	{
	}

	public void run()
	{
		try
		{
			if (writeFn == null)
				writeFn = FileUtils.addBeforeExtention(	spliceFn,
														"_geneNamesAdded");
			if(writeFnNoGeneSymbol==null)
				writeFnNoGeneSymbol=FileUtils.addBeforeExtention(writeFn, "_geneSpliceMerged");

			addGeneNames(	spliceFn,
							new HashMap<String, String[]>(),
							writeFn,
							false);
		} catch (Exception e)
		{
			e.printStackTrace();
		}
	}

	private PerChrIntervalTree<GffElement> getGenome(String gtfFn) throws Exception
	{
		GtfReader gtfReader = new GtfReader(new File(gtfFn));
		PerChrIntervalTree<GffElement> genome = gtfReader.createIntervalTree();
		return genome;
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
			PerChrIntervalTree<GffElement> genome = getGenome(this.gtfFn);
			//add genes to the spliceFile
			//System.out.println("Adding genes to the spliceFile");
			BufferedReader spliceReader = FileUtils.createReader(spliceFN);
			BufferedWriter spliceWriter = FileUtils.createWriter(writeFN);
			BufferedWriter spliceWriter2 = FileUtils.createWriter(writeFnNoGeneSymbol);
			//add the header
			writeFirstLine(spliceReader, spliceFN, afterLine, spliceWriter, spliceWriter2);

			spliceReader.lines().forEach(line ->
			{
				String[] cells = line.split("\t|_");//"_" is for files in which the file junctions are described in 1 column with "_" as separator between aspects of the junction
				String chromLoc = cells[0];// + "\t" + cells[1] + "\t" + cells[2] + "\t" + cells[3] + "\t" + cells[4];
				for (int c = 1; c < 6; c++)
				{
					chromLoc = chromLoc.concat("_");
					chromLoc = chromLoc.concat(cells[c]);
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
								spliceWriter2,
								spliceSiteToGene,
								afterLine);
			});
			spliceWriter2.close();
			spliceReader.close();
			spliceWriter.close();

		} catch (Exception e)
		{
			e.printStackTrace();
		}
		p("GeneNames Merged, file written at:\n" + writeFn);
	};

	private void writeFirstLine(BufferedReader spliceReader, String spliceFN, boolean afterLine, BufferedWriter spliceWriter, BufferedWriter spliceWriter2) throws IOException
	{
		String firstline = spliceReader.readLine();
		if (firstline == null)
		{
			p("File is empty:\t " + spliceFN);
			return;
		}

		if (firstline.contains("First base of the intron (1-based)"))
		{
			if (!afterLine)
			{
				String toWrite = "EnsemblIDs\tGeneSymbols\t" + firstline + "\n";
				spliceWriter.write(toWrite);
				spliceWriter2.write(toWrite);
			}
			else
			{
				String toWrite =firstline + "\t" + "EnsemblIDs\tGeneSymbols" + "\n";
				spliceWriter.write(toWrite);
				spliceWriter2.write(toWrite);
			}
		}
		else if (!firstline.startsWith("\t"))//start the stream at the first line again (in case this is a file that does not have a header)
			spliceReader = FileUtils.createReader(spliceFN);
		else if (!afterLine)
		{
			String toWrite ="EnsemblIds\t" + "GeneSymbols\t" + firstline + "\n";
			spliceWriter.write(toWrite);
			spliceWriter2.write(toWrite);
		}
		else
		{
			String toWrite =firstline + "EnsemblIds\t" + "GeneSymbols\t" + "\n";
			spliceWriter.write(toWrite);
			spliceWriter2.write(toWrite);
		}
	}

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
								PerChrIntervalTree<GffElement> genome,
								BufferedWriter spliceWriter,
								BufferedWriter spliceWriter2,
								HashMap<String, String[]> spliceSiteToGene,
								boolean afterLine)
	{
		try
		{
			//convert junction co-ordinates to exon co-ordinates
			int exonEnd = Integer.parseInt(cells[1])-1;
			int exonStart = Integer.parseInt(cells[2])+1;

			List<GffElement> startOverlappingFeatures = genome.searchPosition(	chromosomeName,
																				exonEnd);
			List<GffElement> endOverlappingFeatures = genome.searchPosition(chromosomeName,
																			exonStart);

			HashMap<String, String[]> genes = new HashMap<String, String[]>();
			genes = addGenesToHashSet(	startOverlappingFeatures,
										genes);
			genes = addGenesToHashSet(	endOverlappingFeatures,
										genes);

			if (genes.size() > 1)
			{
				//if it is overlapping one exon boundary
				HashMap<String, String[]> genesWithMatchingExonBoundaries = new HashMap<String, String[]>();
				genesWithMatchingExonBoundaries = addGeneIfExonBoundaryMatchesToHashSet(startOverlappingFeatures,
																						genesWithMatchingExonBoundaries,
																						exonEnd,
																						exonStart);
				genesWithMatchingExonBoundaries = addGeneIfExonBoundaryMatchesToHashSet(endOverlappingFeatures,
																						genesWithMatchingExonBoundaries,
																						exonEnd,
																						exonStart);
//				System.out.println("genes.size " + genesWithMatchingExonBoundaries.size());
				if(genesWithMatchingExonBoundaries.size()==1)
				{
					String[] genesIDsToAdd = genes.values().iterator().next();
//					System.out.println("gene= " + genesIDsToAdd[0] + "\tSJstart=" + exonEnd);
					writeLineWithIds(	afterLine,
										spliceWriter,
										spliceWriter2,
										genesIDsToAdd,
										line);
				}
			}
			else if (genes.size() == 1)
			{

				String[] genesIDsToAdd = genes.values().iterator().next();
				writeLineWithIds(	afterLine,
									spliceWriter,
									spliceWriter2,
									genesIDsToAdd,
									line);
			}
			else if (genes.size() == 0)
			{
//				addUnnannotated(afterLine,
//								spliceWriter,
//								line,
//								spliceSiteToGene,
//								chromLoc);
			}
		} catch (Exception e)
		{
			e.printStackTrace();
		}
	}

	private HashMap<String, String[]> addGeneIfExonBoundaryMatchesToHashSet(List<GffElement> overlappingFeatures,
																			HashMap<String, String[]> genesWithMatchingExonBoundaries,
																			int exonEnd,
																			int exonStart)
	{
		for (GffElement gffElement : overlappingFeatures)
		{
			if(!gffElement.getFeature().equals("exon"))
				continue;

			//if neither side of the exon is connecting to the junction site
			if(	gffElement.getStart()  !=exonStart&& 
					gffElement.getEnd()  !=exonEnd)
				continue;
			
			String ensemblId = gffElement.getAttributes().get("gene_id");
			String geneSymbol = gffElement.getAttributes().get("gene_name");
			
			if (ensemblId == null)
				ensemblId = "Unannotated";
			if (geneSymbol == null)
				geneSymbol = "Unannotated";
			String[] geneNames = new String[] { ensemblId, geneSymbol };
			genesWithMatchingExonBoundaries.put(	ensemblId,
						geneNames);
		}
		return genesWithMatchingExonBoundaries;
	}

	private HashMap<String, String[]> addGenesToHashSet(List<GffElement> overlappingFeatures,
														HashMap<String, String[]> genes)
	{
		for (GffElement gffElement : overlappingFeatures)
		{
			String ensemblId = gffElement.getAttributes().get("gene_id");
			String geneSymbol = gffElement.getAttributes().get("gene_name");
			if (ensemblId == null)
				ensemblId = "Unannotated";
			if (geneSymbol == null)
				geneSymbol = "Unannotated";
			String[] geneNames = new String[] { ensemblId, geneSymbol };
			if (ensemblId != null)
				;
			genes.put(	ensemblId,
						geneNames);
		}
		return genes;
	}

	private void writeLineWithIds(	boolean afterLine,
									BufferedWriter spliceWriter,
									BufferedWriter spliceWriter2,
									String[] genesIDsToAdd,
									String line) throws IOException
	{
		if (!afterLine)
		{
			spliceWriter.write(genesIDsToAdd[0] + "\t" + genesIDsToAdd[1] + "\t" + line + "\n");
			spliceWriter2.write(genesIDsToAdd[0] + "__" + line + "\n");
		}
		else
			spliceWriter.write(line + "\t" + genesIDsToAdd[0] + "\t" + genesIDsToAdd[1] + "\n");
	}

	private void addUnnannotated(	boolean afterLine,
									BufferedWriter spliceWriter,
									String line,
									HashMap<String, String[]> spliceSiteToGene,
									String chromLoc) throws IOException
	{
		if (!afterLine)
			spliceWriter.write("unAnnotated\tunAnnotated\t" + line + "\n");
		else
			spliceWriter.write(line + "\tunannotated\tunAnnotated" + "\n");
		spliceSiteToGene.put(	chromLoc,
								new String[] { "-", "-" });
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
	
	public String getWriteFn()
	{
		return writeFn;
	}

	public String getGtfFn()
	{
		return gtfFn;
	}

	public void setGtfFn(String gtfFn)
	{
		this.gtfFn = gtfFn;
	}

	public String getWriteFnNoGeneSymbol()
	{
		return writeFnNoGeneSymbol;
	}

	public void setWriteFnNoGeneSymbol(String writeFnNoGeneSymbol)
	{
		this.writeFnNoGeneSymbol = writeFnNoGeneSymbol;
	}
}
