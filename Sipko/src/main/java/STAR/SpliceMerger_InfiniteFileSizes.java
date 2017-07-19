package STAR;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.Array;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Pattern;
import java.util.stream.Collector;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import org.apache.commons.math3.stat.StatUtils;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import MatrixScripts.MatrixString;
import MatrixScripts.MyMatrix;
import Slurm.Slurm;
import Tools.FileUtils;
import Tools.JSONutil;
import Tools.Script;
import Tools.FileSearcher;
import eqtlmappingpipeline.binarymeta.meta.Reader;
import umcg.genetica.containers.Pair;

public class SpliceMerger_InfiniteFileSizes extends Script<SpliceMerger_InfiniteFileSizes> implements Cloneable
{
	// This class does the following:
	// 1. Finds all SJ.out.tab files in a folder and subfolders
	// 2. Combines them into 1 SJ_merged.out.tab file indicating the number of
	// reads overlapping each splice junction in all files
	// 3. Creates additional files with genes that represent these files
	// 4. Adds some other info to the files

	/**
	 * 
	 */
	private static final long serialVersionUID = -7002847670693064746L;
	private String inputFolder_SpliceComment = "/root/directory/,/root/directory2/; MANDATORY // Folder(s) from which all the splice sits should be retrieved (combines all files containing (SJ.out.tab))";
	private String inputFolder_Splice = null;
	private String writeFN_SpliceComment = "/root/directory/mergedSpliceSites.sj.out; MANDATORY // file to which the combined splice junctions are written";
	private String writeFN_Splice = null;// file where all the splice variants merged into 1 file with info on numbe of spliced reads overlapping each in all samples together
	private String writeFN_Splice_2ndPassInputComment = "/root/directory/mergedSpliceSites.sj.out; optional // filenames of the input file for the splice junctions to be used in the second pass. These are all junctions in the fused in the first pass plus those where reads overlapping <= <readCutoff> reads";
	private String writeFN_Splice_2ndPassInput = null;// file where all the splice variants merged into 1 file with info on numbe of spliced reads overlapping each in all samples together
	private String readCutoffComment = "OPTIONAL; // For each splice variant the number of samlpes it is present in is counted. An additional count is made for the number of samples with more then this number of reads overlapping the junction";
	private int readCutoff = 8;
	private String writeFolder_SplicePerGeneComment = "/root/directory/splicePerGene/; MANDATORY // Folder to which splice sites per gene are written";
	private String writeFolder_SplicePerGene = null;// if null becomes: new File(var.writeFN).getParent();
	public String genesToRemove = "ENSG00000244734,ENSG00000188536";//ENSG00000244734,ENSG00000188536 (2 hemoglobine genen)

	// SpliceSitesPerGene variables
	GeneNameAdder geneNameAdder = new GeneNameAdder();
	
	/** Merges the juntion (.sj.out) files into 1 file describing which splice junctions were observed in all .sj.out files in the input folder.
	 * For each splice junction also gives the total number of reads observed in all samples combined. 
	 * 
	 * Paramaters are defined in a .json file
	 */
	
	public void run()
	{
		try
		{
			String sTAR_SpliceFiles = findSpliceFiles();

			Pair<ArrayList<String>, HashMap<String, String[]>> fileNamesCorrected_SpliceSiteToGene = addGeneNamesToSpliceFiles(sTAR_SpliceFiles); // divides reads per splice site by total number of spliced reads for this gene
			
			ArrayList<String> fileNamesCorrected = fileNamesCorrected_SpliceSiteToGene.getLeft();
			HashMap<String, String[]> spliceSiteToGene = fileNamesCorrected_SpliceSiteToGene.getRight();

			//HashTable<SplicesiteName,number of reads overlapping>
			Hashtable<String, SpliceSitesCount> spliceSitesCounts = mergeSpliceFiles(fileNamesCorrected);
			p("Splicejunction reads counted");
			writeSpliceSitesCounts(	spliceSitesCounts,
									spliceSiteToGene);
			p("Splicejunction reads written");

			if (writeFolder_SplicePerGene != null)
			{
				ArrayList<String> geneExpressionFns = new ArrayList<String>();//is filled in writePerGene()
				String writeFolder = writePerGene(	spliceSiteToGene,
													fileNamesCorrected,
													geneExpressionFns);

				p("Reads per gene written to " + writeFolder);
				// mergeSpliceFilesAndZscores(v, fileNamesCorrected, spliceSiteToGene);// writeFN=file with all the splice variants and how often
				writeSpliceSiteToGene(	spliceSiteToGene,
										getWriteFolder_SplicePerGene() + "SpliceSiteToGene.txt"); // the number of times they occur in all the samples together as well as in how many samples they occur

				new FilePerGeneMerger(	geneExpressionFns,
										getWriteFolder_SplicePerGene() + "splicePerGene/expression.txt").run();
			}
		} catch (Exception e)
		{
			e.printStackTrace();
		}
	}

	private String writePerGene(HashMap<String, String[]> spliceSiteToGene,
								ArrayList<String> fileNamesCorrected,
								ArrayList<String> geneExpressionFns) throws FileNotFoundException, IOException
	{
		// get list of genes and the index each splice site will have in the output file
		HashMap<String, ArrayList<String>> geneToSpliceSites = createGeneToSpliceSites(spliceSiteToGene);
		HashMap<String, String> fileNameCorrectedToSampleName = getSampleNames(fileNamesCorrected);

		String folderSplicePerGene = getWriteFolder_SplicePerGene() + "splicePerGene/";
		new File(folderSplicePerGene).mkdir();
		String folderSplicePerGeneExpression = folderSplicePerGene + "Expression/";
		new File(folderSplicePerGeneExpression).mkdir();

		ArrayList<String> keys = new ArrayList<>();
		p("geneToSpliceSites.size = " + geneToSpliceSites.size());

		for (String gene : geneToSpliceSites.keySet())
		{
			keys.add(gene);

		}
		int blockSize = 20000;//write 20.000 genes each round
		for (int b = 0; b <= keys.size() / blockSize; b++)
		{
			int start = b * blockSize;
			int end = (b + 1) * blockSize;
			
			if (end > keys.size())
				end = keys.size();
			// create a writer for each gene and add header to he file
			HashMap<String, BufferedWriter[]> geneWriters = createGeneWriters(	geneToSpliceSites,
																				folderSplicePerGeneExpression,
																				//folderGeneSpliceRatioFN,
																				//geneRatioFns,
																				geneExpressionFns,
																				keys,
																				start,
																				end);

			// load each of the samples into the memory separately

			HashMap<String, HashMap<String, String>> multipletsConvert = new HashMap<String, HashMap<String, String>>();
			int i = 0;
			BufferedWriter writerTest = null;
			for (String fn : fileNamesCorrected)
			{
				if (i % 100 == 0)
					p("File =" + i + "/" + fileNamesCorrected.size());
				HashMap<String, String[]> spliceToData = getSpliceToData(fn);

				// write this sample to the gene files
				int after = 0;
				long file = 0;
				for(String gene: geneWriters.keySet())
				{
					BufferedWriter[] writers = geneWriters.get(gene);
					StringBuilder lineExpression = new StringBuilder();
					lineExpression.append(fileNameCorrectedToSampleName.get(fn));

					for (String splice : geneToSpliceSites.get(gene))
					{
						String[] data = spliceToData.get(splice);
						lineExpression.append("\t");
						if (data == null)
						{
							lineExpression.append("0");
						}
						else
						{

							lineExpression.append(data[6]);// 6 == expression
						}
					}
					lineExpression.append("\n");
					try
					{
						writers[0].write(lineExpression.toString());			
					} catch (Exception e)
					{
						e.printStackTrace();
					}
				}

				i++;
			}
			for(String gene : geneWriters.keySet())
			{
				BufferedWriter[] writers = geneWriters.get(gene);
				writers[0].close();
			}
		}
		return folderSplicePerGene;
	}

	private HashMap<String, BufferedWriter[]> createGeneWriters(HashMap<String, ArrayList<String>> geneToSpliceSites,
																String folderSplicePerGeneExpression,
																ArrayList<String> geneExpressionFns,
																ArrayList<String> keys,
																int start,
																int end) throws FileNotFoundException, IOException
	{
		HashMap<String, BufferedWriter[]> geneWriters = new HashMap<String, BufferedWriter[]>();
		try
		{
			for (int k = start; k < end; k++)
			{
				String gene = keys.get(k);
				if (k % 1000 == 0 || k == keys.size() - 1)
					System.out.println("Gene:" + k + "/" + keys.size());
				BufferedWriter[] writers = new BufferedWriter[2];
				String expressionFn = folderSplicePerGeneExpression + gene + ".txt";
				writers[0] = FileUtils.createWriter(expressionFn);
				geneExpressionFns.add(expressionFn);

				StringBuilder header = new StringBuilder();

				ArrayList<String> spliceSites = geneToSpliceSites.get(gene);
				for (String spliceSite : spliceSites)
				{
					header.append("\t");
					header.append(spliceSite.replace(	"\t",
														"_"));
				}

				header.append("\n");
				writers[0].write(header.toString());
//				writers[1].write(header.toString());

				geneWriters.put(gene,
								writers);
			}
		} catch (Exception e)
		{
			p("Crashed after creating " + geneWriters.size() + " gene writers");
			p("Increased maximum number of files that may be opened simultaneously in the operating system settings");
			p("Exiting");
			e.printStackTrace();
			System.exit(2);
		}
		return geneWriters;
	}

	private HashMap<String, String> getSampleNames(ArrayList<String> fileNamesCorrected)
	{
		HashMap<String, String> fileNameCorrectedToSampleName = new HashMap<String, String>();
		fileNamesCorrected.forEach(fn -> fileNameCorrectedToSampleName.put(	fn,
																			new File(new File(fn).getParent()).getName()));
		return fileNameCorrectedToSampleName;
	}

	private HashMap<String, String[]> getSpliceToData(String fn)
	{
		HashMap<String, String[]> spliceToData = new HashMap<String, String[]>();
		BufferedReader spliceReader;
		try
		{
			spliceReader = FileUtils.createReader(fn);

			spliceReader.lines().forEach(line ->
			{
				String[] eles = line.split("\t");

				String spliceVar = getSpliceVarFromEles(eles);
				spliceToData.put(	spliceVar,
									eles);
				// String spliceVar = spliceBuilder.toString();
				// String annotated = eles[5];
				// int reads = Integer.parseInt(eles[6]);
				// int overhang = Integer.parseInt(eles[8]);
				// int ensemblIDs = eles[9];
				// int geneNames = eles[10];
			});

			spliceReader.close();
		} catch (FileNotFoundException e)
		{
			e.printStackTrace();
		} catch (IOException e)
		{
			e.printStackTrace();
		}
		return spliceToData;
	}

	private String getSpliceVarFromEles(String[] eles)
	{
		StringBuilder spliceBuilder = new StringBuilder();
		spliceBuilder.append(eles[0]);
		spliceBuilder.append("\t");
		spliceBuilder.append(eles[1]);
		spliceBuilder.append("\t");
		spliceBuilder.append(eles[2]);
		spliceBuilder.append("\t");
		spliceBuilder.append(eles[3]);
		spliceBuilder.append("\t");
		spliceBuilder.append(eles[4]);

		return spliceBuilder.toString();
	}

	private HashMap<String, ArrayList<String>> createGeneToSpliceSites(HashMap<String, String[]> spliceSiteToGene)
	{
		HashMap<String, ArrayList<String>> geneToSpliceSites = new HashMap<String, ArrayList<String>>();
		spliceSiteToGene.forEach((	spliceSite,
									gene) ->
		{
			if (gene[0] != null)
			{
				String[] geneNames = gene[0].split(",");

				for (String geneName : geneNames)
				{
					geneName = geneName.toUpperCase();
					ArrayList<String> spliceSites = geneToSpliceSites.get(geneName);
					if (spliceSites == null)
						spliceSites = new ArrayList<String>();
					spliceSites.add(spliceSite);// add the splice site and the index it will have in the output file
					geneToSpliceSites.put(	geneName,
											spliceSites);
				}
			}
		});

		return geneToSpliceSites;
	}

	private void writeSpliceSitesCounts(Hashtable<String, SpliceSitesCount> spliceSitesCounts,
										HashMap<String, String[]> spliceSiteToGene)
	{
		BufferedWriter spliceWriter;
		BufferedWriter spliceWriter2ndPassSites;
		try
		{
			spliceWriter = FileUtils.createWriter(getWriteFN_Splice());
			
			if(getWriteFN_Splice_2ndPassInput()==null)
				setWriteFN_Splice_2ndPassInput(getWriteFN_Splice().replace("SJ_Merged_1stPass.out.tab", "SJ_Merged_2ndPass.input.tab"));//SJ_Merged_1stPass.out.tab
			
			spliceWriter2ndPassSites = FileUtils.createWriter(getWriteFN_Splice_2ndPassInput());//SJ_Merged_1stPass.out.tab
			
			HashMap<String, Integer> removeGenesSpliceCount = new HashMap<String, Integer>();
			
			String spliceWriterGenesRemovedFN = getWriteFN_Splice().replace(".tab",
																			"_" + 0 + "removedGenes.txt");
			if(genesToRemove!=null)
			{
				String[] removeGenes = genesToRemove.split(",");
				
				for (String gene : removeGenes)
					removeGenesSpliceCount.put(	gene,
												0);
				spliceWriterGenesRemovedFN = getWriteFN_Splice().replace(".tab",
																			"_" + removeGenes.length + "removedGenes.txt");
			}

			
			BufferedWriter spliceWriterGenesRemoved = FileUtils.createWriter(spliceWriterGenesRemovedFN);

			String extraWriterFN = getWriteFN_Splice().replace(	".tab",
																"_AdditionalInfo_GeneSymbols.txt");
			BufferedWriter spliceWriter_extra = FileUtils.createWriter(extraWriterFN);

			spliceWriter_extra.write("EnsemblIDs\t" + "GeneSymbols\t" + "chromosome\t" + "First base of the intron (1-based)\t" + "last base of the intron (1-based)\t" + "strand (0: undefined, 1: +, 2: -)\t" + "intron motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5:AT/AC, 6: GT/AT\t" + "0: unannotated, 1: annotated (only if splice junctions database is used)\t" + "number of uniquely mapping reads crossing the junction (in all samples)\t" + "Reserved column\t" + "Maximum overhang detected in all samples\t" + "Number of samples in which this splice junction occurs at least 1 time\t" + "Number of samples in which this splice junction occurs at least " + getReadCutoff() + " times\n");

			for (String spliceVar : spliceSitesCounts.keySet())
			{
				SpliceSitesCount spliceSitesCount = spliceSitesCounts.get(spliceVar);

				StringBuilder line = new StringBuilder();
				line.append(spliceVar);
				line.append("\t");
				line.append(spliceSitesCount.getAnnotated());
				line.append("\t");
				line.append(spliceSitesCount.getReads());
				line.append("\t");
				line.append("0");
				line.append("\t");
				line.append(spliceSitesCount.getMaxOverhang());
				String[] geneIDs = spliceSiteToGene.get(spliceVar);

				String[] ensgGenes = null;
				boolean skip = false;
				if (geneIDs != null && geneIDs[0] != null)
				{
					ensgGenes = geneIDs[0].split(",");

					for (String ensgID : ensgGenes)
					{
						if (removeGenesSpliceCount.containsKey(ensgID))
						{
							removeGenesSpliceCount.put(	ensgID,
														removeGenesSpliceCount.get(ensgID) + 1);
							skip = true;
						}
					}
				}

				if (!skip)
				{
					spliceWriter.write(line.toString() + "\n");
					if(spliceSitesCount.getAnnotated().equals("1") || spliceSitesCount.getReads()>=readCutoff)
						spliceWriter2ndPassSites.write(line.toString() + "\n");
				}
				else
					spliceWriterGenesRemoved.write(line.toString() + "\n");

				StringBuilder line_extra = new StringBuilder();

				if (geneIDs != null)
				{
					line_extra.append(geneIDs[0]);
					line_extra.append("\t");
					line_extra.append(geneIDs[1]);
					line_extra.append("\t");
				}
				else
				{
					line_extra.append("-");
					line_extra.append("\t");
					line_extra.append("-");
					line_extra.append("\t");
				}

				line_extra.append(line);
				line_extra.append("\t");
				line_extra.append(spliceSitesCount.getSamples());
				line_extra.append("\t");
				line_extra.append(spliceSitesCount.getSamplesAboveCutoff());
				line_extra.append("\n");
				spliceWriter_extra.write(line_extra.toString());
			}

			spliceWriter2ndPassSites.close();
			spliceWriter.close();
			spliceWriterGenesRemoved.close();
			spliceWriter_extra.close();
		} catch (FileNotFoundException e)
		{
			e.printStackTrace();
		} catch (IOException e)
		{
			e.printStackTrace();
		}
		p("Merged splicejunction File written to:" + getWriteFN_Splice());
	}

	private Hashtable<String, SpliceSitesCount> mergeSpliceFiles(ArrayList<String> fileNamesCorrected)
	{
		Hashtable<String, SpliceSitesCount> spliceSitesCounts = new Hashtable<String, SpliceSitesCount>();
		for (String fn : fileNamesCorrected)
		{
			spliceSitesCounts = countSpliceSites(	fn,
													spliceSitesCounts);
		}
		return spliceSitesCounts;
	}

	private Hashtable<String, SpliceSitesCount> countSpliceSites(	String fn,
																	Hashtable<String, SpliceSitesCount> spliceSitesCounts)
	{
		try
		{
			BufferedReader spliceReader = FileUtils.createReader(fn);
			spliceReader.lines().forEach(line ->
			{
				String[] eles = line.split("\t");// .*, is any number of characters before a comma; "(?<=)" indicates what pattern has to be present before the first split
				String spliceVar = getSpliceVarFromEles(eles);
				String annotated = eles[5];
				int reads = Integer.parseInt(eles[6]);
				int overhang = Integer.parseInt(eles[8]);
				SpliceSitesCount spliceSitesCount = null;
				if (spliceSitesCounts.containsKey(spliceVar))
					spliceSitesCount = spliceSitesCounts.get(spliceVar);
				else
					spliceSitesCount = new SpliceSitesCount(annotated);
				spliceSitesCount.incrementReads(reads);
				spliceSitesCount.incrementSamples();
				if (reads > getReadCutoff())
					spliceSitesCount.incrementSamplesAboveCutoff();
				if (overhang > spliceSitesCount.getMaxOverhang())
					spliceSitesCount.setMaxOverhang(overhang);

				spliceSitesCounts.put(	spliceVar,
										spliceSitesCount);
			});
			spliceReader.close();
		} catch (FileNotFoundException e)
		{
			e.printStackTrace();
		} catch (IOException e)
		{
			e.printStackTrace();
		}

		return spliceSitesCounts;
	}

	private Pair<ArrayList<String>, HashMap<String, String[]>> addGeneNamesToSpliceFiles(String sTAR_SpliceFiles) // divides reads per splice site by total number of
																												// spliced reads for this gene
																												throws FileNotFoundException, IOException
	{
		BufferedReader readerFiles = FileUtils.createReader(sTAR_SpliceFiles);
		// create two hashes to return results in
		HashMap<String, String[]> spliceSiteToGene = new HashMap<String, String[]>();
		ArrayList<String> fileNamesWithGeneNames = new ArrayList<String>();
		// Add gene names to each file
		readerFiles.lines().forEach(spliceFile -> fileNamesWithGeneNames.add(addGeneNamesOneFile(spliceFile,
																										spliceSiteToGene)));// normalize(file,spliceAvgStdev,zScoreFolder,spliceHash)
		readerFiles.close();
		// return output
		Pair<ArrayList<String>, HashMap<String, String[]>> pair = new Pair<ArrayList<String>, HashMap<String, String[]>>(	fileNamesWithGeneNames,
																															spliceSiteToGene);
		p("SpliceSites detected:\t" + spliceSiteToGene.size());
		
		return pair;
	}

	//Add gene names to each file
	private String addGeneNamesOneFile(	String spliceFile,
													HashMap<String, String[]> spliceSiteToGene)
	{
		String outputFN = null;
		try
		{
			outputFN = geneNameAdder.run(	spliceFile,
													spliceSiteToGene);
		} catch (FileNotFoundException e)
		{
			e.printStackTrace();
		} catch (IOException e)
		{
			e.printStackTrace();
		}
		return outputFN;
	}

	private String findSpliceFiles() throws Exception
	{
		String writeFN_Filenames = getWriteFolder_SplicePerGene() + "SJ_FileNames.txt";
		FileSearcher spliceFinder = new FileSearcher();

		if (this.jsonFN != null)
			spliceFinder.setJsonFN(FileUtils.removeExtention(writeFN_Filenames) + "FileSearcher.json");
		spliceFinder.setFolders(getInputFolder_Splice());
		spliceFinder.setWriteName(writeFN_Filenames);
		spliceFinder.setSearchStrings(new String[] { "SJ.out.tab" });
		spliceFinder.setForbiddenStrings(new String[] { "_GeneSymbols" });
		spliceFinder.run();
		return writeFN_Filenames;
	}

	private void writeSpliceSiteToGene(	HashMap<String, String[]> startEndToGene,
										String writeFN)
	{
		try
		{
			BufferedWriter writer = FileUtils.createWriter(writeFN);
			startEndToGene.forEach((spliceVar,
									ensg_GeneSymbol) ->
			{
				String line = spliceVar.replace("\t",
												"_")
						+ "\t" + ensg_GeneSymbol[0] + "\t" + ensg_GeneSymbol[1];
				try
				{
					writer.write(line + "\n");
				} catch (Exception e)
				{
					e.printStackTrace();
				}
			});

			writer.close();
		} catch (Exception e)
		{
			e.printStackTrace();
		}
	}

	void mergeSpliceSitesPerGene(	STAR_Pipeline v,
									HashMap<String, String[]> startEndToGene,
									String folderPerSplice,
									ArrayList<String> sTAR_SpliceFiles)
	{
		HashMap<String, ArrayList<String>> geneToSpliceFiles = createGeneToSpliceFileshash(	v,
																							startEndToGene,
																							folderPerSplice);

		String folderSplicePerGene = getWriteFolder_SplicePerGene() + "splicePerGene/";
		new File(folderSplicePerGene).mkdir();
		String folderSplicePerGeneExpression = folderSplicePerGene + "spliceExpression/";
		new File(folderSplicePerGeneExpression).mkdir();
		

		String[] rowNames = new String[sTAR_SpliceFiles.size()];
		sTAR_SpliceFiles.toArray(rowNames);
		String[] rowNames2 = Stream.of(rowNames).map(rowName -> new File(new File(rowName).getParent()).getName()).collect(Collectors.toList()).stream().toArray(size -> new String[size]);// how																																														// on
																																															// earth
																																															// do
																																															// I
																																															// return
																																															// an
																																															// array
																																															// of
																																															// strings
																																															// more
																																															// smoothly
																																															// from
																																															// a
																																															// stream?
		geneToSpliceFiles.forEach((	gene,
									files) -> writeGeneSplicevariantsToFile(gene,
																			files,
																			folderSplicePerGeneExpression,
																			rowNames2,
																			startEndToGene));
	}

	private void writeGeneSplicevariantsToFile(	String gene,
												ArrayList<String> files,
												String folderSplicePerGeneExpression,
												String[] rowNames,
												HashMap<String, String[]> startEndToGene)
	{

		int n = 0;
		for (String file : files)
			if (new File(file).exists())
				n++;

		MatrixString mergeMatrix = new MatrixString(rowNames.length,
													n);
		mergeMatrix.set0();
		mergeMatrix.rowNames = rowNames;;
		Hashtable<String, Integer> rowHash = mergeMatrix.getRowHash();
		int outCol = 0;
		for (int f = 0; f < files.size(); f++)
		{
			String spliceFN = files.get(f);
			if (!new File(spliceFN).exists())
				continue;
			MatrixString spliceFile = new MatrixString(spliceFN);

			String spliceVarName = new File(spliceFN).getName().replace(".txt",
																		"");
			String spliceVar = spliceVarName.replace(	"_multiGene",
														"").replace("_",
																	"\t");
			mergeMatrix.setColHeader(	outCol,
										spliceVarName);

			String[] geneIDs = startEndToGene.get(spliceVar);
			int pos = 0;
			if (geneIDs != null)
			{
				String[] geneName = geneIDs[0].split(",");
				pos = Arrays.asList(geneName).indexOf(gene);
			}

			for (int r = 0; r < spliceFile.rows(); r++)
			{
				int outRow = rowHash.get(spliceFile.rowNames[r]);
				mergeMatrix.values[outRow][outCol] = spliceFile.values[r][0];;
			}
			outCol++;
		}
		String writeGeneSpliceReadsFN = folderSplicePerGeneExpression + gene + ".txt";

		mergeMatrix.write(writeGeneSpliceReadsFN);
	}

	private HashMap<String, ArrayList<String>> createGeneToSpliceFileshash(	STAR_Pipeline v,
																			HashMap<String, String[]> startEndToGene,
																			String folderPerSplice)
	{
		HashMap<String, ArrayList<String>> geneToSpliceFiles = new HashMap<String, ArrayList<String>>();
		startEndToGene.forEach((location,
								names) ->
		{
			ArrayList<String> spliceFiles = null;
			if (names[0] != null)
			{
				String[] geneIdentifiers = names[0].split(",");// gene symbol
																// (names[0] =
																// ensembl id)
				for (int g = 0; g < geneIdentifiers.length; g++)
				{
					String geneIdentifier = geneIdentifiers[g];
					if (geneToSpliceFiles.containsKey(geneIdentifier))
					{
						spliceFiles = geneToSpliceFiles.get(geneIdentifier);
					}
					else
						spliceFiles = new ArrayList<String>();
					String spliceFN = folderPerSplice + location.replace(	"\t",
																			"_")
							+ ".txt";
					if (geneIdentifiers.length > 1)// this splice variant
													// overlaps multiple genes
						spliceFN = spliceFN.replace(".txt",
													"_multiGene.txt");
					spliceFiles.add(spliceFN);
					geneToSpliceFiles.put(	geneIdentifier,
											spliceFiles);
				}
			}
		});
		return geneToSpliceFiles;
	}

	public String getInputFolder_Splice()
	{
		return inputFolder_Splice;
	}

	public void setInputFolder_Splice(String inputFolder_Splice)
	{
		this.inputFolder_Splice = inputFolder_Splice;
	}

	public String getWriteFN_Splice()
	{
		return writeFN_Splice;
	}

	public void setWriteFN_Splice(String writeFN_Splice)
	{
		this.writeFN_Splice = writeFN_Splice;
	}

	public int getReadCutoff()
	{
		return readCutoff;
	}

	public void setReadCutoff(int readCutoff)
	{
		this.readCutoff = readCutoff;
	}

	public String getWriteFolder_SplicePerGene()
	{
		return writeFolder_SplicePerGene;
	}

	public void setWriteFolder_SplicePerGene(String writeFolder_Splice)
	{
		this.writeFolder_SplicePerGene = writeFolder_Splice;
	}

	public String getWriteFN_Splice_2ndPassInput()
	{
		return writeFN_Splice_2ndPassInput;
	}

	public void setWriteFN_Splice_2ndPassInput(String writeFN_Splice_2ndPassInput)
	{
		this.writeFN_Splice_2ndPassInput = writeFN_Splice_2ndPassInput;
	}
}
