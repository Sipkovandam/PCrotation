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

import Kallisto.Slurm;
import PCA.Matrix;
import PCA.MatrixString;
import Tools.FileUtils;
import Tools.JSONutil;
import Tools.SearchFilesInDirectories;
import eqtlmappingpipeline.binarymeta.meta.Reader;
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

	static void run(Variables v) throws IOException 
	{
		run(v, true);
	}
	
	static void run(Variables v, boolean writePerGene) throws IOException {
		String sTAR_SpliceFiles = findSpliceFiles(v);

		Pair<ArrayList<String>, HashMap<String, String[]>> fileNamesCorrected_SpliceSiteToGene = correctForGeneExpression(v, sTAR_SpliceFiles); //divides reads per splice site by total number of spliced reads for this gene
		ArrayList<String> fileNamesCorrected = fileNamesCorrected_SpliceSiteToGene.getLeft();
		HashMap<String, String[]> spliceSiteToGene = fileNamesCorrected_SpliceSiteToGene.getRight();
		
		Hashtable<String, SpliceSitesCount> spliceSitesCounts = mergeSpliceFiles(v, fileNamesCorrected);
		System.out.println("Splicejunction reads counted");
		writeSpliceSitesCounts(v,spliceSitesCounts, spliceSiteToGene);
		System.out.println("Splicejunction reads written");
		
		if(writePerGene)
		{
			String writeFolder = writePerGene(v, spliceSiteToGene,fileNamesCorrected);
			
			System.out.println("Reads per gene written to " + writeFolder);
			//mergeSpliceFilesAndZscores(v, fileNamesCorrected, spliceSiteToGene);// writeFN=file with  all the splice variants and how often
			writeSpliceSiteToGene(spliceSiteToGene,v.getWriteFolder_Splice()+"SpliceSiteToGene.txt");													  // they occur in all the samples together as well as in how manysamples they occur
		}
	}

	private static String writePerGene(Variables v , HashMap<String, String[]> spliceSiteToGene, ArrayList<String> fileNamesCorrected) throws FileNotFoundException, IOException 
	{
		//get list of genes and the index each splice site will have in the output file
		HashMap<String, ArrayList<String>> geneToSpliceSites = createGeneToSpliceSites(spliceSiteToGene);
		
		HashMap<String, String> fileNameCorrectedToSampleName = getSampleNames(fileNamesCorrected);

		String folderSplicePerGene=  v.getWriteFolder_Splice() + "splicePerGene/";
		new File(folderSplicePerGene).mkdir();
		String folderSplicePerGeneExpression = folderSplicePerGene+"spliceExpression/";
		new File(folderSplicePerGeneExpression).mkdir();
		String folderGeneSpliceRatioFN = folderSplicePerGene+"spliceRatio/";
		new File(folderGeneSpliceRatioFN).mkdir();
		
		//create a writer for each gene and add header to he file
		HashMap<String, BufferedWriter[]> geneWriters = createGeneWriters(geneToSpliceSites, folderSplicePerGeneExpression, folderGeneSpliceRatioFN);
		
		//load each of the samples into the memory separately
		
		HashMap<String,HashMap<String,String>> multipletsConvert = new HashMap<String,HashMap<String,String>>();
		for(String fn : fileNamesCorrected)
		{
			HashMap<String, String[]> spliceToData = getSpliceToData(fn);
			
			//write this sample to the gene files
			geneWriters.forEach((gene,writers) ->
			{
				StringBuilder lineRatio = new StringBuilder();
				StringBuilder lineExpression = new StringBuilder(); 
				lineRatio.append(fileNameCorrectedToSampleName.get(fn));
				lineExpression.append(fileNameCorrectedToSampleName.get(fn));
				
				for(String splice : geneToSpliceSites.get(gene))
				{
					String[] data = spliceToData.get(splice);
					lineRatio.append("\t");
					lineExpression.append("\t");
					if(data == null)
					{
						lineRatio.append("0");
						lineExpression.append("0");
					}
					else
					{
						String ratio = getRightRatio(data, multipletsConvert, gene);
							
						lineRatio.append(ratio);//11 == ratio
						lineExpression.append(data[6]);//6 == expression
					}
				}
				lineRatio.append("\n");
				lineExpression.append("\n");
				try {
					writers[0].write(lineExpression.toString());
					writers[1].write(lineRatio.toString());
				} catch (Exception e) {	e.printStackTrace();}
			});
			
			
		}
		for(BufferedWriter[] writers : geneWriters.values())
		{
			writers[0].close();
			writers[1].close();
		}
		return folderSplicePerGene;
	}

	private static String getRightRatio(String[] data, HashMap<String, HashMap<String, String>> multipletsConvert, String gene) 
	{
		String ratio = data[11];
		if(ratio.contains(","))
		{
			String geneNames = data[10];
			HashMap<String,String> geneToRatio = multipletsConvert.get(geneNames);
			if(geneToRatio == null)
			{
				String[] ratios = ratio.split(",");
				String[] gNs = data[10].split(",");
				geneToRatio = new HashMap<String,String>();
				for(int g = 0; g < gNs.length; g++)
				{
					geneToRatio.put(gNs[g], ratios[g]);
				}
				multipletsConvert.put(geneNames, geneToRatio);
			}
			ratio = geneToRatio.get(gene);
		}
		return ratio;
	}

	private static HashMap<String, BufferedWriter[]> createGeneWriters(HashMap<String, ArrayList<String>> geneToSpliceSites, String folderSplicePerGeneExpression, String folderGeneSpliceRatioFN) throws FileNotFoundException, IOException {
		HashMap<String, BufferedWriter[]> geneWriters = new HashMap<String, BufferedWriter[]>();
		for(String gene : geneToSpliceSites.keySet())
		{
			BufferedWriter[] writers = new BufferedWriter[2];
			writers[0] = FileUtils.createWriter(folderSplicePerGeneExpression+gene+".txt");
			writers[1] = FileUtils.createWriter(folderGeneSpliceRatioFN+gene+".txt");

			StringBuilder header = new StringBuilder();
			
			ArrayList<String> spliceSites = geneToSpliceSites.get(gene);
			for(String spliceSite : spliceSites)
			{
				header.append("\t");
				header.append(spliceSite.replace("\t", "_"));
			}

			header.append("\n");
			writers[0].write(header.toString());
			writers[1].write(header.toString());
			
			geneWriters.put(gene, writers);
		}
		return geneWriters;
}
private static HashMap<String, String> getSampleNames(ArrayList<String> fileNamesCorrected) 
	{
		HashMap<String, String> fileNameCorrectedToSampleName = new HashMap<String, String>();
		fileNamesCorrected.forEach(fn -> fileNameCorrectedToSampleName.put(fn, new File(new File(fn).getParent()).getName()));
		return fileNameCorrectedToSampleName;
	}
private static HashMap<String, String[]> getSpliceToData(String fn) 
	{
		HashMap<String, String[]> spliceToData = new HashMap<String, String[]>();
		BufferedReader spliceReader;
		try {
			spliceReader = FileUtils.createReader(fn);
		
			spliceReader.lines().forEach(line -> 
			{
				String[] eles = line.split("\t");
				
				String spliceVar = getSpliceVarFromEles(eles);
				spliceToData.put(spliceVar, eles);
//				String spliceVar =  spliceBuilder.toString();
//				String annotated = eles[5];
//				int reads = Integer.parseInt(eles[6]);
//				int overhang = Integer.parseInt(eles[8]);
//				int ensemblIDs = eles[9];
//				int geneNames = eles[10];
//				double ratios = double.parseDouble(eles[11]);
				
				
			});
			
			spliceReader.close();
		} catch (FileNotFoundException e) {e.printStackTrace();} catch (IOException e) {e.printStackTrace();}
		return spliceToData;
	}

	private static String getSpliceVarFromEles(String[] eles) {
		StringBuilder spliceBuilder = new StringBuilder();
		spliceBuilder.append(eles[0]);spliceBuilder.append("\t");
		spliceBuilder.append(eles[1]);spliceBuilder.append("\t");
		spliceBuilder.append(eles[2]);spliceBuilder.append("\t");
		spliceBuilder.append(eles[3]);spliceBuilder.append("\t");
		spliceBuilder.append(eles[4]);
		
		return spliceBuilder.toString();
}

	private static HashMap<String, ArrayList<String>> createGeneToSpliceSites(HashMap<String, String[]> spliceSiteToGene) 
	{
		HashMap<String, ArrayList<String>> geneToSpliceSites = new HashMap<String, ArrayList<String>>();
		spliceSiteToGene.forEach((spliceSite,gene)->
		{
			if(gene[1] !=null)
			{
				String[] geneSymbols = gene[1].split(",");
				
				for(String geneSymbol : geneSymbols)
				{
					ArrayList<String> spliceSites = geneToSpliceSites.get(geneSymbol);
					if(spliceSites == null)
						spliceSites = new ArrayList<String>();
					spliceSites.add(spliceSite);//add the splice site and the index it will have in the output file
					geneToSpliceSites.put(geneSymbol, spliceSites);
				}
			}
		});
		
		return geneToSpliceSites;
	}

	private static void writeSpliceSitesCounts(Variables v, Hashtable<String, SpliceSitesCount> spliceSitesCounts, HashMap<String, String[]> spliceSiteToGene) {
		BufferedWriter spliceWriter;
		try {
			spliceWriter = FileUtils.createWriter(v.getWriteFN_Splice());
		
		String extraWriterFN = v.getWriteFN_Splice().replace(".tab", "_AdditionalInfo_GeneSymbols.txt");
		BufferedWriter spliceWriter_extra = FileUtils.createWriter(extraWriterFN);
		spliceWriter_extra.write("EnsemblIDs\t" + "GeneSymbols\t" + "chromosome\t"
				+ "First base of the intron (1-based)\t" + "last base of the intron (1-based)\t"
				+ "strand (0: undefined, 1: +, 2: -)\t"
				+ "intron motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5:AT/AC, 6: GT/AT\t"
				+ "0: unannotated, 1: annotated (only if splice junctions database is used)\t"
				+ "number of uniquely mapping reads crossing the junction (in all samples)\t"
				+ "Reserved column\t"
				+ "Maximum overhang detected in all samples\t"
				+ "Number of samples in which this splice junction occurs at least 1 time\t"
				+ "Number of samples in which this splice junction occurs at least " + v.getReadCutoff() + " times\n");

		spliceSitesCounts.forEach((spliceVar,spliceSitesCount) -> 
		{
			StringBuilder line = new StringBuilder();
			line.append(spliceVar);line.append("\t");
			line.append(spliceSitesCount.getAnnotated());line.append("\t");
			line.append(spliceSitesCount.getReads());line.append("\t");
			line.append("0");line.append("\t");
			line.append(spliceSitesCount.getMaxOverhang());
			
			try {
				spliceWriter.write(line.toString()+"\n");
			} catch (Exception e) {e.printStackTrace();}
			
			StringBuilder line_extra = new StringBuilder();
			String[] geneIDs = spliceSiteToGene.get(spliceVar);
			if(geneIDs!=null)
			{
				line_extra.append(geneIDs[0]);line_extra.append("\t");
				line_extra.append(geneIDs[1]);line_extra.append("\t");
			}
			else
			{
				line_extra.append("-");line_extra.append("\t");
				line_extra.append("-");line_extra.append("\t");
			}
				
			line_extra.append(line);line_extra.append("\t");
			line_extra.append(spliceSitesCount.getSamples());line_extra.append("\t");
			line_extra.append(spliceSitesCount.getSamplesAboveCutoff());line_extra.append("\n");
			
			try {
				spliceWriter_extra.write(line_extra.toString());
			} catch (Exception e) {	e.printStackTrace();}
		});

		spliceWriter.close();
		spliceWriter_extra.close();
		} catch (FileNotFoundException e) {e.printStackTrace();} catch (IOException e) {e.printStackTrace();}
		System.out.println("Merged splicejunction File written to:" + v.getWriteFN_Splice());
	}

	private static Hashtable<String, SpliceSitesCount> mergeSpliceFiles(Variables v, ArrayList<String> fileNamesCorrected) 
	{
		Hashtable<String, SpliceSitesCount> spliceSitesCounts = new Hashtable<String, SpliceSitesCount>();
		for(String fn : fileNamesCorrected)
		{
			spliceSitesCounts = countSpliceSites(v, fn, spliceSitesCounts);
		}
		return spliceSitesCounts;
	}

	private static Hashtable<String, SpliceSitesCount> countSpliceSites(Variables v, String fn,
			Hashtable<String, SpliceSitesCount> spliceSitesCounts) 
	{
		try {
			BufferedReader spliceReader = FileUtils.createReader(fn);
			spliceReader.lines().forEach(line -> 
			{				
				String[] eles = line.split("\t");//.*, is any number of characters before a comma; "(?<=)" indicates what pattern has to be present before the first split
				String spliceVar = getSpliceVarFromEles(eles);
				String annotated = eles[5];
				int reads = Integer.parseInt(eles[6]);
				int overhang = Integer.parseInt(eles[8]);
				SpliceSitesCount spliceSitesCount = null;
				if(spliceSitesCounts.containsKey(spliceVar))
					spliceSitesCount = spliceSitesCounts.get(spliceVar);
				else
					spliceSitesCount = new SpliceSitesCount(annotated);
				spliceSitesCount.incrementReads(reads);
				spliceSitesCount.incrementSamples();
				if(reads> v.getReadCutoff())
					spliceSitesCount.incrementSamplesAboveCutoff();
				if(overhang>spliceSitesCount.getMaxOverhang())
					spliceSitesCount.setMaxOverhang(overhang);
				
				spliceSitesCounts.put(spliceVar, spliceSitesCount);
			});
			spliceReader.close();
		} catch (FileNotFoundException e) {	e.printStackTrace();} catch (IOException e) 
		{e.printStackTrace();}
		
		return spliceSitesCounts;
	}

	private static Pair<ArrayList<String>, HashMap<String, String[]>> correctForGeneExpression(Variables v,//divides reads per splice site by total number of spliced reads for this gene
			String sTAR_SpliceFiles) throws FileNotFoundException, IOException {
		BufferedReader readerFiles = FileUtils.createReader(sTAR_SpliceFiles);
		HashMap<String, String[]> spliceSiteToGene = new HashMap<String, String[]>();
		ArrayList<String> fileNamesNormalized = new ArrayList<String>();
		readerFiles.lines().forEach(spliceFile -> fileNamesNormalized.add(normalize(v, spliceFile, spliceSiteToGene)));// normalize(file,spliceAvgStdev,zScoreFolder,spliceHash)
		readerFiles.close();
		Pair<ArrayList<String>, HashMap<String, String[]>> pair = new Pair<ArrayList<String>, HashMap<String, String[]>>(
				fileNamesNormalized, spliceSiteToGene);
		return pair;
	}

	private static String normalize(Variables v, String spliceFile, HashMap<String, String[]> spliceSiteToGene) {
		String outputFN = null;
		try {
			outputFN = NormalizeForGeneExpression.run(v, spliceFile, spliceSiteToGene);
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

	private static void writeSpliceSiteToGene(HashMap<String, String[]> startEndToGene, String writeFN) {
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
