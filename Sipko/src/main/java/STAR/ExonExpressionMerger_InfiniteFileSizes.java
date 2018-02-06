package STAR;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Stream;

import Tools.FileSearcher;
import Tools.FileUtils;
import Tools.Script;
import umcg.genetica.containers.Pair;

public class ExonExpressionMerger_InfiniteFileSizes extends Script<ExonExpressionMerger_InfiniteFileSizes>
{
	/**
	 * Also corrects exon expression for exon length when calculating exon expression ratios (relative to whole gene expression)
	 */
	private static final long serialVersionUID = -8940680788150581125L;
	private String searchFolder = null;//semitransient
	private String writeFolder = null;//semitransient
	private String writeMergedFn = null;//semitransient
	private String ensgToGeneSymbolFn = null;//semitransient

	public void run()
	{
		log(searchFolder);
		try
		{
			init();
			// find all the exon files (featureCounts output files)
			String featureCountFilesFn = findExonFiles();
			// read geneID conversion table
			HashMap<String, String> ensgToGeneSymbol = FileUtils.readStringStringHash(ensgToGeneSymbolFn);

			// create 1 writer for each gene and get all exonsIDs per gene
			String writeFolderExpression = getWriteFolder() + "exonPerGene/";

			ArrayList<String> geneExpressionFns = new ArrayList<String>();//is filled in createWritersPerGene()

			List<Object> geneToExon_exonToLength = getGeneToExon_getExonToLength(	featureCountFilesFn,
																ensgToGeneSymbol,
																writeFolderExpression,
																geneExpressionFns);
			//			HashMap<String, BufferedWriter> writerPerGene = (HashMap<String, BufferedWriter>) returnObjects.get(0);
			HashMap<String, ArrayList<String>> geneToExons = (HashMap<String, ArrayList<String>>) geneToExon_exonToLength.get(0);
			HashMap<String, Double> exonToLength = (HashMap<String, Double>) geneToExon_exonToLength.get(1);

			new File(writeFolderExpression).mkdir();

			int blockSize = 20000;//write 20.000 genes each round
			
			String[] genes = new String[geneToExons.size()];
			int i = 0;
			for(String gene : geneToExons.keySet())
			{
				genes[i]=gene;
				i++;
			}
			for (int b = 0; b <= geneToExons.size() / blockSize; b++)
			{
				int start = b * blockSize;
				int end = (b + 1) * blockSize;

				if (end > geneToExons.size())
					end = geneToExons.size();

				HashMap<String, BufferedWriter> writerPerGene = new HashMap<String, BufferedWriter>();
				
				for(int g = start; g < end; g++)
				{
					String gene = genes[g];
					writerPerGene.put(	gene,
										writeHeader(gene,
										            geneToExons.get(gene),
													writeFolderExpression,
													geneExpressionFns,
													start,
													end));
				}

				// go through all the featureCounts files one by one and write the
				// results per sample to results per gene files
				log("Creating expression files");
				BufferedReader fcFnReader = FileUtils.createReader(featureCountFilesFn);
				AtomicInteger counter = new AtomicInteger(0);
				fcFnReader.lines().forEach(fn -> writeToGeneFiles(	fn,
																	writerPerGene,
																	geneToExons,
																	ensgToGeneSymbol,
																	counter,
																	geneExpressionFns,
																	writeFolderExpression));
				log(counter.get() + " files merged");
				// close all writers
				for (BufferedWriter writer : writerPerGene.values())
					writer.close();

				// go through all the expression files (so each gene)
				log("Creating ratio files");

				String folderName = getWriteFolder() + "Ratios/";
				new File(folderName).mkdir();

				ArrayList<String> geneRatioFns = new ArrayList<String>();

				File[] expressionFiles = new File(writeFolderExpression).listFiles();
				for (File expressionFile : expressionFiles)
				{
					BufferedReader reader = FileUtils.createReader(expressionFile.getAbsolutePath());
					String expressionFN = folderName + expressionFile.getName();
					BufferedWriter writer = FileUtils.createWriter(expressionFN);
					geneRatioFns.add(expressionFN);
					// copy the first line
					String headerLine = reader.readLine();
					String[] exons = headerLine.split("\t");
					writer.write(headerLine + "\n");
					reader.lines().forEach(line -> calcAndWriteRatios(	line,
																		writer,
																		exonToLength,
																		exons));
					writer.close();
				}
				new FilePerGeneMerger(	geneRatioFns,
				                      	writeMergedFn).run();
			}
		} catch (Exception e)
		{
			e.printStackTrace();
		}

	}

	private void calcAndWriteRatios(String line,
									BufferedWriter writer,
									HashMap<String, Double> exonToLength,
									String[] exons)
	{
		try
		{
			String[] valuesStr = line.split("\t");
			double[] values = new double[valuesStr.length];
			double sum = 0;

			for (int v = 1; v < valuesStr.length; v++)
			{
				values[v] = Double.parseDouble(valuesStr[v]) / exonToLength.get(exons[v]);
				sum += values[v];
			}

			StringBuilder writeLine = new StringBuilder();
			writeLine.append(valuesStr[0]);

			for (int v = 1; v < values.length; v++)
			{
				writeLine.append("\t");
				double ratio = values[v] / sum;
				if (Double.isNaN(ratio))
					ratio = -1;
				writeLine.append(ratio);
			}
			writeLine.append("\n");
			writer.write(writeLine.toString());
		} catch (Exception e)
		{
			System.out.println("line =\t" + line);
			String[] valuesStr = line.split("\t");
			double[] values = new double[valuesStr.length];
			double sum = 0;
			for (int v = 1; v < valuesStr.length; v++)
			{
				System.out.println("exons[v]=\t" + exons[v]);
				System.out.println("exonToLength.get(exons[v] =\t" + exonToLength.get(exons[v]) + "\t" + exonToLength.size());
				System.out.println("exonToLength.get(exons[v] =\t" + exonToLength.get(exonToLength.keySet().iterator().next()));
				values[v] = Double.parseDouble(valuesStr[v]) / exonToLength.get(exons[v]);
				sum += values[v];
				System.out.println("lvaluesStr[v]" + valuesStr[v]);
				System.out.println("exonToLength.get(exons[v]) =\t" + exonToLength.get(exons[v]));
			}
			
			e.printStackTrace();
		}
	}

	private void init()
	{
		if (writeFolder != null)
		{
			File writeFolderFile = new File(writeFolder);
			if (!writeFolderFile.exists())
			{
				File parent = writeFolderFile.getParentFile();
				if (parent.exists())
					writeFolderFile.mkdir();
				else
				{
					log("Parent folder does not exist: \t" + parent.getAbsolutePath());
					log("Exiting");
					System.exit(2);

				}
			}
		}

		if(writeMergedFn==null)
			writeMergedFn=writeFolder + "expression.txt.gz";
	}

	private void writeToGeneFiles(	String fn,
									HashMap<String, BufferedWriter> writers,
									HashMap<String, ArrayList<String>> geneToExons,
									HashMap<String, String> ensgToGeneSymbol,
									AtomicInteger i,
									ArrayList<String> geneExpressionFns,
									String writeFolderExpression)
	{

		if (i.getAndIncrement() % 100 == 0)
			log(i.get() + " files merged");
		BufferedReader exonFileReader;
		HashMap<String, String> exonToData = new HashMap<String, String>();
		try
		{
			exonFileReader = FileUtils.createReader(fn);
			//get the data (reads per exon) for each exon
			exonFileReader.lines().skip(2).forEach(exonLine -> addDataToExon(	exonLine,
																				writers,
																				ensgToGeneSymbol,
																				exonToData));

			// writePerGene
			writers.forEach((	gene,
								writer) ->
			{
				String sample = new File(new File(fn).getParent()).getName();
				StringBuilder lineBuilder = new StringBuilder();
				lineBuilder.append(sample);
				ArrayList<String> exons = geneToExons.get(gene);
				for (String exon : exons)
				{
					String reads = exonToData.get(exon);
					lineBuilder.append("\t");
					lineBuilder.append(reads);
				}
				try
				{
					lineBuilder.append("\n");
					writer.write(lineBuilder.toString());
				} catch (Exception e)
				{
					e.printStackTrace();
				}
			});
		} catch (FileNotFoundException e)
		{
			e.printStackTrace();
		} catch (IOException e)
		{
			e.printStackTrace();
		}
	}

	private void addDataToExon(	String exonLine,
								HashMap<String, BufferedWriter> writers,
								HashMap<String, String> ensgToGeneSymbol,
								HashMap<String, String> exonToData)
	{
		String[] eles = exonLine.split("\t");
		String ensgName = eles[0];
		String reads = eles[6];
		String exon = getExonName(eles);

		exonToData.put(	exon,
						reads);

	}

	private List<Object> getGeneToExon_getExonToLength(	String featureCountFilesFn,
												HashMap<String, String> ensgToGeneSymbol,
												String writeFolderExpression,
												ArrayList<String> geneExpressionFns) throws FileNotFoundException, IOException
	{
		String firstFn = FileUtils.createReader(featureCountFilesFn).readLine();
		BufferedReader exonFileReader = FileUtils.createReader(firstFn);
		HashMap<String, ArrayList<String>> geneToExons = new HashMap<>();
		HashMap<String, Double> exonToLength = new HashMap<>();
		exonFileReader.lines().skip(2).forEach(line -> fillExonHashMaps(line,
																		geneToExons,
																		ensgToGeneSymbol,
																		exonToLength));

		return r(	geneToExons,
					exonToLength);
	}

	private BufferedWriter writeHeader(	String gene,
										ArrayList<String> exons,
										String writeFolderExpression,
										ArrayList<String> geneExpressionFns,
										int start,
										int end)
	{
		BufferedWriter writer = null;
		try
		{
			String fn = writeFolderExpression + gene + ".txt";
			writer = FileUtils.createWriter(fn);
			geneExpressionFns.add(fn);
			for (int e = 0; e < exons.size(); e++)
				writer.write("\t" + exons.get(e));
			writer.write("\n");
		} catch (FileNotFoundException e)
		{
			e.printStackTrace();
		} catch (IOException e)
		{
			e.printStackTrace();
		}
		return writer;
	}

	private void fillExonHashMaps(	String line,
									HashMap<String, ArrayList<String>> exonsPerGene,
									HashMap<String, String> ensgToGeneSymbol,
									HashMap<String, Double> exonToLength)
	{
		try
		{
			String[] eles = line.split("\t");

			String exon = getExonName(eles);
			exonToLength.put(	exon,
								Double.parseDouble(eles[5]));
			String ensgName = eles[0];

			String geneSymbol = ensgToGeneSymbol.get(ensgName);
			
			if (geneSymbol == null)
				return;
			String geneID = geneSymbol.toUpperCase();
			// geneID=ensgName;

			ArrayList<String> exons = exonsPerGene.get(geneID);
			if (exons == null)
				exons = new ArrayList<String>();
			exons.add(exon);
			exonsPerGene.put(	geneID,
								exons);

		} catch (Exception e)
		{
			String[] eles = line.split("\t");
			String ensgName = eles[0];
			System.out.println("ensgName=\t" + ensgName);
			String geneSymbol = ensgToGeneSymbol.get(ensgName);
			System.out.println("ensgname=\t" + geneSymbol);
			e.printStackTrace();
		}

	}

	private String getExonName(String[] eles)
	{
		String chr = eles[1];
		String start = eles[2];
		String end = eles[3];
		StringBuilder exonBuilder = new StringBuilder();
		exonBuilder.append(chr);
		exonBuilder.append("_");
		exonBuilder.append(start);
		exonBuilder.append("_");
		exonBuilder.append(end);
		return exonBuilder.toString();
	}

	private String findExonFiles()
	{
		String writeFN_Filenames = getWriteFolder() + "featureCounts_FileNames.txt";
		FileSearcher spliceFinder = new FileSearcher();

		if (this.jsonFN != null)
			spliceFinder.setJsonFN(FileUtils.removeExtention(writeFN_Filenames) + "FileSearcher.config");
		spliceFinder.setFolders(getSearchFolder());
		spliceFinder.setWriteName(writeFN_Filenames);
		spliceFinder.setSearchStrings(new String[] { "featureCounts.out" });
		spliceFinder.setForbiddenStrings((new String[] { "featureCounts.out.summary" }));
		spliceFinder.run();
		return writeFN_Filenames;
	}

	public String getWriteFolder()
	{
		return FileUtils.makeFolderNameEndWithSlash(writeFolder);
	}

	public void setWriteFolder(String writeFolder)
	{
		this.writeFolder = writeFolder;
	}

	public String getSearchFolder()
	{
		return searchFolder;
	}

	public void setSearchFolder(String searchFolder)
	{
		this.searchFolder = searchFolder;
	}

	public String getEnsgToGeneSymbolFn()
	{
		return ensgToGeneSymbolFn;
	}

	public void setEnsgToGeneSymbolFn(String ensgToGeneSymbolFn)
	{
		this.ensgToGeneSymbolFn = ensgToGeneSymbolFn;
	}

	public String getWriteMergedFn()
	{
		return writeMergedFn;
	}

	public void setWriteMergedFn(String writeMergedFn)
	{
		this.writeMergedFn = writeMergedFn;
	}
}
