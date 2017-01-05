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
import java.util.stream.Collector;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import org.apache.commons.math3.stat.StatUtils;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import Kallisto.Slurm;
import PCA.Matrix;
import PCA.MatrixString;
import Tools.FileUtils;
import Tools.JSONutil;
import Tools.SearchFilesInDirectories;
import umcg.genetica.containers.Pair;

public class SpliceSites {
	// This class does the following:
	// 1. Finds all SJ.out.tab files in a folder and subfolders
	// 2. Combines them into 1 SJ_merged.out.tab file indicating the number of
	// reads overlapping each splice junction in all files
	// 3. Creates additional files with genes that represent these files
	// 4. Adds z-scores and some other stats to the files

	public static void main(String args[]) throws Exception {
		// args = new
		// String[]{"json=E:/Groningen/Test/STAR/STAR/SpliceSites.json"};
		Variables v = checkArgs(args);
		v.writeVars();

		run(v);
	}

	static void run(Variables v) throws IOException {
		String sTAR_SpliceFiles = findSpliceFiles(v);

		Pair<ArrayList<String>, HashMap<String, String[]>> fileNamesCorrected_startEndToGene = correctForTotalGeneExpression(v, sTAR_SpliceFiles); //divides reads per splice site by total number of spliced reads for this gene
		ArrayList<String> fileNamesCorrected = fileNamesCorrected_startEndToGene.getLeft();
		HashMap<String, String[]> startEndToGene = fileNamesCorrected_startEndToGene.getRight();
		mergeSpliceFilesAndZscores(v, fileNamesCorrected, startEndToGene);// writeFN=file with  all the splice variants and how often
																		  // they occur in all the samples together as well as in how manysamples they occur
	}

	private static Pair<ArrayList<String>, HashMap<String, String[]>> correctForTotalGeneExpression(Variables v,//divides reads per splice site by total number of spliced reads for this gene
			String sTAR_SpliceFiles) throws FileNotFoundException, IOException {
		BufferedReader readerFiles = FileUtils.createReader(sTAR_SpliceFiles);
		HashMap<String, String[]> startEndToGene = new HashMap<String, String[]>();
		ArrayList<String> fileNamesNormalized = new ArrayList<String>();
		readerFiles.lines().forEach(spliceFile -> fileNamesNormalized.add(normalize(v, spliceFile, startEndToGene)));// normalize(file,spliceAvgStdev,zScoreFolder,spliceHash)
		readerFiles.close();
		Pair<ArrayList<String>, HashMap<String, String[]>> pair = new Pair<ArrayList<String>, HashMap<String, String[]>>(
				fileNamesNormalized, startEndToGene);
		return pair;
	}

	private static String normalize(Variables v, String spliceFile, HashMap<String, String[]> startEndToGene) {
		String outputFN = null;
		try {
			outputFN = NormalizeForGeneExpression.run(v, spliceFile, startEndToGene);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return outputFN;
	}

	private static String findSpliceFiles(Variables v) throws IOException {
		String writeFN_Filenames = v.getWriteFolder_Splice() + "SJ_FileNames.txt";
		System.out.println("v.getWriteFolder_Splice() = " + v.getWriteFolder_Splice());
		System.out.println("v.getInputFolder_Splice() = " + v.getInputFolder_Splice());
		SearchFilesInDirectories.main(new String[] { "foldername=" + v.getInputFolder_Splice(),
				"writefn=" + writeFN_Filenames, "searchstrings=" + "SJ.out.tab" });
		return writeFN_Filenames;
	}

	private static void mergeSpliceFilesAndZscores(Variables v, ArrayList<String> sTAR_SpliceFiles,
			HashMap<String, String[]> startEndToGene) throws FileNotFoundException, IOException {
		System.out.println("Summing splice reads in all files to one file");
		HashMap<String, int[]> spliceHash = new HashMap<String, int[]>();
		HashMap<String, ArrayList<SpliceCount>> splicesPerSampleHash = new HashMap<String, ArrayList<SpliceCount>>();// has the fragments per sample (in SpliceCountobjects) for splice variant (keys)
		AtomicInteger nRefSamples = new AtomicInteger();
		sTAR_SpliceFiles.forEach(file -> nRefSamples.set(addToHash(v, file, spliceHash, splicesPerSampleHash,
				v.getExcludeFromReferenceSamplesContaining().split(","), nRefSamples.get())));
		// write the hash
		BufferedWriter spliceWriter = FileUtils.createWriter(v.getWriteFN_Splice());
		String extraWriterFN = v.getWriteFN_Splice().replace(".tab", "_AdditionalInfo_GeneSymbols.txt");
		BufferedWriter spliceWriter_extra = FileUtils.createWriter(extraWriterFN);
		spliceWriter_extra.write("EnsemblIDs\t" + "GeneSymbols\t" + "chromosome\t"
				+ "First base of the intron (1-based)\t" + "last base of the intron (1-based)\t"
				+ "strand (0: undened, 1: +, 2: -)\t"
				+ "intron motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5:AT/AC, 6: GT/AT\t"
				+ "0: unannotated, 1: annotated (only if splice junctions database is used)\t"
				+ "number of uniquely mapping reads crossing the junction (in all samples)\t"
				+ "number of multi-mapping reads crossing the junction\t" + "maximum spliced alignment overhang\t"
				+ "Number of samples in which this splice junction occurs at least 1 time\t"
				+ "Number of samples in which this splice junction occurs at least " + v.getReadCutoff() + " times\t"
				+ "Average (all samples): % of spliced reads mapping to this splice variant relative to all those mapping to all spliced reads overlapping this gene\t"
				+ "Standard deviation (all samples): % of spliced reads mapping to this splice variant relative to all those mapping to all spliced reads overlapping this gene\n");

		System.out.println(
				"calculate the averages and standard deviations of the number of reads overlapping each transcript in each sample");
		HashMap<String, double[]> spliceAvgStdev = new HashMap<String, double[]>();
		splicesPerSampleHash.forEach((transcript, countsPerSample) -> spliceAvgStdev.put(transcript,
				getAvgStdev(countsPerSample, nRefSamples.get(), transcript)));

		HashMap<String, int[]> spliceHashPresentInMost = new HashMap<String, int[]>();

		String folderPerSplice = v.getWriteFolder_Splice() + "expressionPerSplice/";
		new File(folderPerSplice).mkdir();

		// nodig:
		// sample+splice --> rank
		// hashtable<SampleName,SpliceVariant> --> rank
		HashMap<String, HashMap<String, Integer>> perSamplePerSpliceToRank = new HashMap<String, HashMap<String, Integer>>(
				sTAR_SpliceFiles.size());
		System.out.println("Writing All Splice Sites");
		Runtime runtime = Runtime.getRuntime();
		NumberFormat format = NumberFormat.getInstance();

		StringBuilder sb = new StringBuilder();
		long maxMemory = runtime.maxMemory();
		long allocatedMemory = runtime.totalMemory();
		long freeMemory = runtime.freeMemory();

		sb.append("free memory: " + format.format(freeMemory / 1024 / 1024) + "\n");
		sb.append("allocated memory: " + format.format(allocatedMemory / 1024 / 1024) + "\n");
		sb.append("max memory: " + format.format(maxMemory / 1024 / 1024) + "\n");
		sb.append("total free memory: " + format.format((freeMemory + (maxMemory - allocatedMemory)) / 1024 / 1024)
				+ "\n");
		System.out.println(sb.toString());

		AtomicInteger counter = new AtomicInteger();
		spliceHash.forEach((transcript, spliceCounts) -> {
			boolean presentInMost = writeAllSplicesSites(v, transcript, spliceCounts, spliceWriter, spliceWriter_extra,
					splicesPerSampleHash, spliceAvgStdev, sTAR_SpliceFiles.size(), folderPerSplice,
					perSamplePerSpliceToRank, spliceHash.size(), counter, startEndToGene);

			if (presentInMost)
				spliceHashPresentInMost.put(transcript, spliceCounts);
		});

		maxMemory = runtime.maxMemory();
		allocatedMemory = runtime.totalMemory();
		freeMemory = runtime.freeMemory();
		sb.setLength(0);
		sb.append("free memory: " + format.format(freeMemory / 1024 / 1024) + "\n");
		sb.append("allocated memory: " + format.format(allocatedMemory / 1024 / 1024) + "\n");
		sb.append("max memory: " + format.format(maxMemory / 1024 / 1024) + "\n");
		sb.append("total free memory: " + format.format((freeMemory + (maxMemory - allocatedMemory)) / 1024 / 1024)
				+ "\n");
		System.out.println(sb.toString());

		spliceWriter.close();
		spliceWriter_extra.close();

		// System.out.println("Adding GeneNames And Counting Reads Per Splice
		// variant");
		// SpliceSitesPerGene.addGeneNamesAndCountReadsPerSplice(v,
		// extraWriterFN,
		// null,
		// new HashMap<String, String[]>(),
		// null,
		// false);

		// add z-scores to file
		String zScoreFolder = v.getWriteFolder_Splice() + "zScoresPerSample/";
		new File(zScoreFolder).mkdir();
//		System.out.println("Adding ZscoresPerSample");
//		sTAR_SpliceFiles.forEach(file -> writeZscoresPerSample(v, file, spliceAvgStdev, zScoreFolder, spliceHash,
//				spliceHashPresentInMost, startEndToGene, perSamplePerSpliceToRank));

		mergeSpliceSitesPerGene(v, startEndToGene, folderPerSplice, sTAR_SpliceFiles);
		
		writeStartEndToGene(startEndToGene,v.getWriteFolder_Splice()+"StartEndToGene.txt");
	}

	private static void writeStartEndToGene(HashMap<String, String[]> startEndToGene, String writeFN) {
		try {
			BufferedWriter writer = FileUtils.createWriter(writeFN);
			startEndToGene.forEach((spliceVar,ensg_GeneSymbol)->
			{
				String line = spliceVar.replace("\t", "_")+"\t"+ensg_GeneSymbol[0]+"\t"+ensg_GeneSymbol[1];
				try {
					writer.write(line+"\n");
				} catch (Exception e) {e.printStackTrace();}
			});
			
			
			writer.close();
		} catch (Exception e) {	e.printStackTrace();} 	
	}

	private static double[] getAvgStdev(ArrayList<SpliceCount> countsPerSample, int nSamples, String transcript) {
		double[] values = new double[nSamples];
		int v = 0;
		for (int i = 0; i < countsPerSample.size(); i++) {
			if (!countsPerSample.get(i).getIncludeInReference())
				continue;
			values[v] = countsPerSample.get(i).getRelativeAbundances()[0];
			v++;
		}

		double average = StatUtils.mean(values);
		double stdev = Math.sqrt(StatUtils.variance(values));
		double[] stats = new double[] { average, stdev };

		return stats;
	}

	private static boolean writeAllSplicesSites(Variables v, String spliceVar, int[] spliceCounts,
			BufferedWriter spliceWriter, BufferedWriter spliceWriter_extra,
			HashMap<String, ArrayList<SpliceCount>> perSampelHash, HashMap<String, double[]> spliceAvgStdev,
			int nSamples, String folderPerSplice, HashMap<String, HashMap<String, Integer>> perSamplePerSpliceToRank,
			int nSpliceSites, AtomicInteger counter, HashMap<String, String[]> startEndToGene) {
		StringBuilder lineBuilder = new StringBuilder();
		lineBuilder.append(spliceVar);	lineBuilder.append("\t");
		lineBuilder.append(spliceCounts[0]);	lineBuilder.append("\t");
		lineBuilder.append(spliceCounts[1]);	lineBuilder.append("\t");
		lineBuilder.append(spliceCounts[2]);
		String line = lineBuilder.toString();
		String chromLoc = spliceVar;
		counter.set(counter.get() + 1);
		if (counter.get() % 1000 == 0)
			System.out.println("SpliceVar:" + counter.get() + "/" + nSpliceSites);
		try {
			StringBuilder lineBuilder_extra = new StringBuilder();

			spliceWriter.write(line + "\n");
			String ensemblIDs=startEndToGene.get(chromLoc)[0];
			String geneSymbol_IDs=startEndToGene.get(chromLoc)[1];
			lineBuilder_extra.append(ensemblIDs);	lineBuilder_extra.append("\t");
			lineBuilder_extra.append(geneSymbol_IDs);	lineBuilder_extra.append("\t");
			lineBuilder_extra.append(line);	lineBuilder_extra.append("\t");
			lineBuilder_extra.append(spliceCounts[3]);	lineBuilder_extra.append("\t");
			lineBuilder_extra.append(spliceCounts[4]);	lineBuilder_extra.append("\t");
			double[] avgStdev = spliceAvgStdev.get(spliceVar);
			lineBuilder_extra.append(avgStdev[0]);	lineBuilder_extra.append("\t");
			lineBuilder_extra.append(avgStdev[1]);

			spliceWriter_extra.write(lineBuilder_extra.toString() + "\n");
			if (perSampelHash.get(spliceVar).size() >= nSamples * v.getMinPercentageSamplesPresentForSpliceToBeWrittenToOwnFile())
			{
				String writeFN = folderPerSplice + spliceVar.replace("\t", "_") + ".txt";
				if(geneSymbol_IDs!=null && geneSymbol_IDs.contains(","))//this splice variant overlaps multiple genes
				{
					writeFN=writeFN.replace(".txt", "_multiGene.txt");
				}
				writeSpliceVarFile(v, perSampelHash.get(spliceVar), spliceVar,
						writeFN, perSamplePerSpliceToRank,
						nSpliceSites);
			}
			if (spliceCounts[3] >= v.getMinPercentage() * nSamples)
				return true;
		} catch (IOException e) {
			e.printStackTrace();
		}
		return false;
	}

	private static int[] writeSpliceVarFile(Variables v, ArrayList<SpliceCount> spliceCountsPerSample, String spliceVar,
			String spliceVariantWriterFN, HashMap<String, HashMap<String, Integer>> perSamplePerSpliceToRank,
			int nSpliceSites) {
		int[] outlierNumbers = new int[2];

		try (FileWriter spliceVariantWriter = new FileWriter(new File(spliceVariantWriterFN))) {// FileUtils.createWriter(spliceVariantWriterFN);){
			spliceVariantWriter.write("Sample\tReadsMappedToThisSpliceJunction\tRelativeNumberOfReadsComparedToReadsSplicedToSameGene\n");

			Collections.sort(spliceCountsPerSample, (SpliceCount s1, SpliceCount s2) -> Double
					.compare(s1.getRelativeAbundances()[0], s2.getRelativeAbundances()[0]));//this is used to more quickly identify the ranks of genes below

			int smaller = 0;
			
			
			for (int s = 0; s < spliceCountsPerSample.size(); s++)// for each sample
			{
				SpliceCount spliceCount = spliceCountsPerSample.get(s);
				
				double[] abundances = spliceCount.getRelativeAbundances();
				String abuncancesStr="";
				abuncancesStr += abundances[0];
				for(double abundance : spliceCount.getRelativeAbundances())
					abuncancesStr+= ","+abundance;
				String line = spliceCount.getSampleName() + "\t" + spliceCount.getReadsOverlapping() + "\t"
						+ abuncancesStr;
				try {
					spliceVariantWriter.write(line + "\n");
				} catch (Exception e) {
					e.printStackTrace();
				}

				if (spliceCount.getReadsOverlapping() > v.getReadCutoff() && spliceCount.getIncludeInReference())
					smaller++;

				// put the rank of this sample for this variant based on its
				// relative abundance (compared to other splice variants of the
				// same gene) into a hash
				if (perSamplePerSpliceToRank.containsKey(spliceCount.getSampleName())) {
					HashMap<String, Integer> spliceToRank = perSamplePerSpliceToRank.get(spliceCount.getSampleName());
					spliceToRank.put(spliceVar, smaller);
					perSamplePerSpliceToRank.put(spliceCount.getSampleName(), spliceToRank);
				} else {
					HashMap<String, Integer> spliceToRank = new HashMap<String, Integer>(nSpliceSites, 0.5f);
					spliceToRank.put(spliceVar, smaller);
					perSamplePerSpliceToRank.put(spliceCount.getSampleName(), spliceToRank);
				}
			}
			//
			// spliceCountsPerSample.stream().forEach(spliceCount -> //for each
			// sample
			// {
			// String line =
			// spliceCount.getSampleName()+"\t"+spliceCount.getReadsOverlapping()+"\t"+spliceCount.getRelativeAbundance();
			// try {
			// spliceVariantWriter.write(line+"\n");
			// } catch (Exception e) {e.printStackTrace();}

			// if(spliceCount.getReadsOverlapping()>v.getReadCutoff())
			// smaller.getAndIncrement();
			//
			// //put the rank of this sample for this variant based on its
			// relative abundance (compared to other splice variants of the same
			// gene) into a hash
			// if(perSamplePerSpliceToRank.containsKey(spliceCount.getSampleName()))
			// {
			// HashMap<String,Integer> spliceToRank=
			// perSamplePerSpliceToRank.get(spliceCount.getSampleName());
			// spliceToRank.put(spliceVar, smaller.get());
			// perSamplePerSpliceToRank.put(spliceCount.getSampleName(),
			// spliceToRank);
			// }
			// else
			// {
			// HashMap<String,Integer> spliceToRank= new
			// HashMap<String,Integer>();
			// spliceToRank.put(spliceVar, smaller.get());
			// perSamplePerSpliceToRank.put(spliceCount.getSampleName(),
			// spliceToRank);
			// }
			// });

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return outlierNumbers;
	}

	private static void writeZscoresPerSample(Variables v, String filename, HashMap<String, double[]> spliceAvgStdev,
			String zScoreFolder, HashMap<String, int[]> spliceHash, HashMap<String, int[]> spliceHashPresentInMost,
			HashMap<String, String[]> startEndToGene,
			HashMap<String, HashMap<String, Integer>> perSamplePerSpliceToRank) {
		BufferedReader readerSpliceFile;
		BufferedWriter writerSpliceFile;
		BufferedWriter writerSpliceFile2;
		try {
			readerSpliceFile = FileUtils.createReader(filename);
			writerSpliceFile = FileUtils.createWriter(filename.replace(".tab", "_zScoresAdded.out.tab"));
			String writeName2 = new File(new File(filename).getParent()).getName() + "SJ_zScoresAdded.out.tab";
			writerSpliceFile2 = FileUtils.createWriter(zScoreFolder + writeName2);
			String header = "chromosome\t" + "First base of the intron (1-based)\t"
					+ "last base of the intron (1-based)\t" + "strand (0: undened, 1: +, 2: -)\t"
					+ "intron motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5:AT/AC, 6: GT/AT\t"
					+ "0: unannotated, 1: annotated (only if splice junctions database is used)\t"
					+ "number of uniquely mapping reads crossing the junction\t"
					+ "number of multi-mapping reads crossing the junction\t" + "maximum spliced alignment overhang\t"
					+ "Ensembl Gene IDs\t" + "Gene Symbol\t" + "% of spliced reads mapping to this gene\t"
					+ "Number of samples in which this splice junction occurs at least 1 time\t"
					+ "Number of REFRENCE samples in which this splice junction occurs at least " + v.getReadCutoff()
					+ " times\t" + "zScores\t" + "Average In All REFRENCE Samples\t" + "Stdev In All REFRENCE Samples\t"
					+ "Number of samples (including this sample) with a relatively lower abundance for this splice variant, compared to other samples (Only samples with > "
					+ v.getReadCutoff() + " are considered)";
			writerSpliceFile.write(header + "\n");
			writerSpliceFile2.write(header + "\n");

			HashSet<String> inFile = new HashSet<String>();
			String sampleName = new File(filename).getParentFile().getName();
			readerSpliceFile.lines().forEach(line -> {
				String[] eles = line.split("\t");
				StringBuilder spliceVarBuilder = new StringBuilder();
				spliceVarBuilder.append(eles[0]);
				for (int i = 1; i < 6; i++) {
					spliceVarBuilder.append("\t");
					spliceVarBuilder.append(eles[i]);
				}
				String spliceVar = spliceVarBuilder.toString();
				inFile.add(spliceVar);

				double[] avgStdev = spliceAvgStdev.get(spliceVar);
				double average = avgStdev[0], stdev = avgStdev[1];
				String[] ratios = eles[11].split(",");
				double zScore = (Double.parseDouble(ratios[0]) - average) / stdev;
				String newLine = line + "\t" + spliceHash.get(spliceVar)[3] + "\t" + spliceHash.get(spliceVar)[4] + "\t"
						+ zScore + "\t" + average + "\t" + stdev + "\t"
						+ perSamplePerSpliceToRank.get(sampleName).get(spliceVar);
				try {
					writerSpliceFile.write(newLine + "\n");
					writerSpliceFile2.write(newLine + "\n");
				} catch (Exception e) {
					e.printStackTrace();
				}
			});

			// System.out.println("infile = " + inFile.size());
			// System.out.println(spliceHashPresentInMost.size());
			spliceHashPresentInMost.forEach((spliceVar, spliceCounts) -> {
				if (!inFile.contains(spliceVar)) {
					double[] avgStdev = spliceAvgStdev.get(spliceVar);
					double average = avgStdev[0], stdev = avgStdev[1];
					double zScore = (0 - average) / stdev;
					String chromLoc = spliceVar;
					String newLine = spliceVar + "\t" + 0 + "\t" + 0 + "\t" + spliceHash.get(spliceVar)[2] + "\t"
							+ startEndToGene.get(chromLoc)[0] + "\t" + startEndToGene.get(chromLoc)[1] + "\t" + 0 + "\t"
							+ spliceHash.get(spliceVar)[3] + "\t" + spliceHash.get(spliceVar)[4] + "\t" + zScore + "\t"
							+ average + "\t" + stdev+ "\t"+ perSamplePerSpliceToRank.get(sampleName).get(spliceVar);
					try {
						writerSpliceFile.write(newLine + "\n");
						writerSpliceFile2.write(newLine + "\n");
					} catch (Exception e) {
						e.printStackTrace();
					}
				}
			});

			writerSpliceFile.close();
			writerSpliceFile2.close();
			readerSpliceFile.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static int addToHash(Variables v,
								 String filename, 
								 HashMap<String, int[]> spliceHash,
								 HashMap<String, ArrayList<SpliceCount>> perSampleHash, 
								 String[] excludeSampleStrings, int nRefSamples) {
		BufferedReader readerSpliceFile;
		try {
			readerSpliceFile = FileUtils.createReader(filename);
			String sampleName = new File(filename).getParentFile().getName();
			System.out.println("s="+ nRefSamples + " SampleName = " + sampleName);
			boolean refSample = Stream.of(excludeSampleStrings).filter(ex -> filename.contains(ex))
					.collect(Collectors.toList()).size() == 0;
			if (refSample == true)
				nRefSamples++;
			readerSpliceFile.lines()
					.forEach(line -> addToSpliceHash(v, line, spliceHash, sampleName, perSampleHash, refSample));

			readerSpliceFile.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return nRefSamples;
	}

	private static void addToSpliceHash(Variables v, String line, HashMap<String, int[]> spliceCountHash,
			String sampleName, HashMap<String, ArrayList<SpliceCount>> perSpliceSampleHash, boolean isRef) {
		String[] eles = line.split("\t");
		StringBuilder spliceVarBuilder = new StringBuilder();
		spliceVarBuilder.append(eles[0]);spliceVarBuilder.append("\t");
		spliceVarBuilder.append(eles[1]);spliceVarBuilder.append("\t");
		spliceVarBuilder.append(eles[2]);spliceVarBuilder.append("\t");
		spliceVarBuilder.append(eles[3]);spliceVarBuilder.append("\t");
		spliceVarBuilder.append(eles[4]);spliceVarBuilder.append("\t");
		spliceVarBuilder.append(eles[5]);
		String spliceVar = spliceVarBuilder.toString();
		int observations = Integer.parseInt(eles[6]);
		int maxOverhang = Integer.parseInt(eles[8]);
		
		String[] abundancesStr= eles[11].split(",");
		double[] relativeAbundances = new double[abundancesStr.length];
		for(int i = 0; i < relativeAbundances.length; i++)
			relativeAbundances[i] = Double.parseDouble(abundancesStr[i]);
		// System.out.println(eles[1]);
		// if(eles[1].contains("3980067") && eles[2].contains("3980513"))
		// System.out.println("Observations: " + eles[6]);

		if (!spliceCountHash.containsKey(spliceVar)) {
			int[] values = new int[5];
			values[0] = observations;// how many reads overlap this junction
			values[1] = 0;// how often it is observed in multimapped reads
							// (multimapped reads are ignored using the settings
							// I use)
			values[2] = maxOverhang;// max overhang
			values[3] = 1;// number of samples it is observed in
			if (observations > v.getReadCutoff() && isRef)// values[4] will hold the number of samples with 8 (2^3) or more reads overlapping this junction
				values[4] = 1;
			spliceCountHash.put(spliceVar, values);

			// put the samplename and the number of reads that map to it into
			// this hash (these will be the ending columns of the file)
			// String samples =
			// sampleName+"\t"+observations+"\t"+relativeAbundance;

			SpliceCount spliceCount = new SpliceCount(sampleName, observations, relativeAbundances, isRef);
			ArrayList<SpliceCount> spliceCounts = new ArrayList<SpliceCount>();
			spliceCounts.add(spliceCount);
			perSpliceSampleHash.put(spliceVar, spliceCounts);
		} else {

			int[] values = spliceCountHash.get(spliceVar);
			values[0] += observations;// how many reads overlap this junction
			values[1] += 0;// how often it is observed in multimapped reads
							// (multimapped reads are ignored using the settings
							// I use)
			if (maxOverhang > values[2])
				values[2] = maxOverhang;// max overhang
			values[3]++;// number of samples it is observed in
			if (observations > v.getReadCutoff() && isRef)
				values[4]++;
			spliceCountHash.put(spliceVar, values);

			// put the samplename and the number of reads that map to it into
			// this hash (these will be the ending columns of the file)
			// String samples = perSpliceSampleHash.get(spliceVar);
			ArrayList<SpliceCount> spliceCounts = perSpliceSampleHash.get(spliceVar);
			SpliceCount spliceCount = new SpliceCount(sampleName, observations, relativeAbundances, isRef);
			spliceCounts.add(spliceCount);
			perSpliceSampleHash.put(spliceVar, spliceCounts);
		}
	}

	static Variables checkArgs(String[] args) {
		if (System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko") && args.length < 1)
			return new Variables();
		if (args.length < 1) {
			System.out.println("Script requires the following argumetns:\n" + "1. ...=<...txt> - ...\n");
			System.exit(1);
		}

		for (int a = 0; a < args.length; a++) {
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase()) {
			case "json":
				return Variables.readVars(value);
			default:
				System.out.println("Incorrect argument supplied:\n" + args[a] + "\nexiting");
				System.exit(1);
			}
		}
		return null;
	}

	static void mergeSpliceSitesPerGene(Variables v, HashMap<String, String[]> startEndToGene, String folderPerSplice,
			ArrayList<String> sTAR_SpliceFiles) {
		HashMap<String, ArrayList<String>> geneToSpliceFiles = createGeneToSpliceFileshash(v, startEndToGene,
				folderPerSplice);

		String folderSplicePerGene=  v.getWriteFolder_Splice() + "splicePerGene/";
		new File(folderSplicePerGene).mkdir();
		String folderSplicePerGeneExpression = folderSplicePerGene+"spliceExpression/";
		new File(folderSplicePerGeneExpression).mkdir();
		String folderGeneSpliceRatioFN = folderSplicePerGene+"spliceRatio/";
		new File(folderGeneSpliceRatioFN).mkdir();
		
		String[] rowNames = new String[sTAR_SpliceFiles.size()];
		sTAR_SpliceFiles.toArray(rowNames);
		String[] rowNames2 = Stream.of(rowNames).map(rowName -> new File(new File(rowName).getParent()).getName())
				.collect(Collectors.toList()).stream().toArray(size -> new String[size]);// how on earth do I return an array of strings more smoothly from a stream?
		geneToSpliceFiles
				.forEach((gene, files) -> writeGeneSplicevariantsToFile(gene, files, folderSplicePerGeneExpression, folderGeneSpliceRatioFN, rowNames2, startEndToGene));
	}

	private static void writeGeneSplicevariantsToFile(String gene, 
													  ArrayList<String> files, 
													  String folderSplicePerGeneExpression,
													  String folderGeneSpliceRatioFN,
													  String[] rowNames, 
													  HashMap<String, String[]> startEndToGene) 
	{
		
		int n = 0;
		for(String file : files)
			if(new File(file).exists())
				n++;
		
		MatrixString mergeMatrix = new MatrixString(rowNames.length, n);
		mergeMatrix.set0();
		MatrixString mergeMatrixRatios = new MatrixString(rowNames.length, n);
		mergeMatrixRatios.set0();
		mergeMatrix.rowNames = rowNames;
		mergeMatrixRatios.rowNames = rowNames;
		Hashtable<String, Integer> rowHash = mergeMatrix.getRowHash();
		int outCol = 0;
		for (int f = 0; f < files.size(); f++) {
			String spliceFN = files.get(f);
			if(!new File(spliceFN).exists())
				continue;
			MatrixString spliceFile = new MatrixString(spliceFN);
			
			String spliceVarName = new File(spliceFN).getName().replace(".txt", "");
			String spliceVar=spliceVarName.replace("_multiGene", "").replace("_", "\t");
			mergeMatrix.setColHeader(outCol, spliceVarName);
			mergeMatrixRatios.setColHeader(outCol, spliceVarName);
			
			String[] geneIDs = startEndToGene.get(spliceVar);
			int pos = 0;
			if(geneIDs!=null)
			{
				String[] geneSymbol= geneIDs[1].split(",");
				pos = Arrays.asList(geneSymbol).indexOf(gene);
			}
			
			for (int r = 0; r < spliceFile.rows(); r++) {
				int outRow = rowHash.get(spliceFile.rowNames[r]);
				mergeMatrix.values[outRow][outCol] = spliceFile.values[r][0];
				mergeMatrixRatios.values[outRow][outCol] = spliceFile.values[r][1].split(",")[pos];
			}
			outCol++;
		}
		String writeGeneSpliceReadsFN = folderSplicePerGeneExpression+ gene + ".txt";
		
		mergeMatrix.write(writeGeneSpliceReadsFN);
		String writeGeneSpliceRatioFN = folderGeneSpliceRatioFN+ gene + ".txt";
		mergeMatrixRatios.write(writeGeneSpliceRatioFN);
		//System.out.println("File written to: " + writeGeneSpliceRatioFN);
	}

	private static HashMap<String, ArrayList<String>> createGeneToSpliceFileshash(Variables v,
			HashMap<String, String[]> startEndToGene, String folderPerSplice) {
		HashMap<String, ArrayList<String>> geneToSpliceFiles = new HashMap<String, ArrayList<String>>();
		startEndToGene.forEach((location, names) -> {
			ArrayList<String> spliceFiles = null;
			if(names[1] != null)
			{
				String[] geneIdentifiers = names[1].split(",");// gene symbol (names[0] = ensembl id)
				for(int g = 0; g < geneIdentifiers.length; g++)
				{
					String geneIdentifier = geneIdentifiers[g];
					if (geneToSpliceFiles.containsKey(geneIdentifier)) {
						spliceFiles = geneToSpliceFiles.get(geneIdentifier);
					} 
					else
						spliceFiles = new ArrayList<String>();
					String spliceFN = folderPerSplice + location.replace("\t", "_") + ".txt"; 
					if(geneIdentifiers.length>1)//this splice variant overlaps multiple genes
						spliceFN=spliceFN.replace(".txt", "_multiGene.txt");
					spliceFiles.add(spliceFN);
					geneToSpliceFiles.put(geneIdentifier, spliceFiles);
				}
			}
		});
		return geneToSpliceFiles;
	}
}
