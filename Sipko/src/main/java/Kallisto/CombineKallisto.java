package Kallisto;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import JuhaPCA.PCA;
import MatrixScripts.MyMatrix;
import Slurm.SumTranscriptsToGenes;
import Tools.FileUtils;
import Tools.JSONutil;
import Tools.Script;
import Tools.FileSearcher;

public class CombineKallisto extends Script<CombineKallisto>{

	String kallistoOutputFolder = null;
	String writeFolder = null;
	int kallistoColumn = 2;// 2 is read counts, 3 is TPM values in
									// kallisto output
	String transcriptsToGenesFN = null;
	double threshold = 0.7;
	int minPercentageFeaturesExpressed = 0;// features are either
													// transcripts or genes
	String mappingPercentagesFn = null;// "E:/Groningen/Test/JSON/ServerTest/Kallisto/mappingPerSample.txt";
	String tsvThatshouldBeThereFN = null;// "E:/Groningen/Test/JSON/ServerTest/Kallisto/scriptNumberToFiles.txt";//should have the .tsv filenames in the first column (including path)
	String genesFN = null;//optional; writeFn for counts per gene
	
	@Override
	public void run() 
	{
		try
		{
			if (writeFolder == null)
				writeFolder = new File(kallistoOutputFolder).getAbsolutePath() + "/";
			
			if(mappingPercentagesFn == null)
				mappingPercentagesFn= FileUtils.makeFolderNameEndWithSlash(writeFolder)+"mappingPercentages.txt";
			
			String combindedFN = writeFolder + "counts.txt.gz";
			if(genesFN == null)
			{
				genesFN=FileUtils.replaceEnd(combindedFN, "_GENES.txt.gz");
			}
			
			this.writeConfig();
			
			if (tsvThatshouldBeThereFN != null) {
				ArrayList<String> tsvThatshouldBeThere = FileUtils.readArrayList(tsvThatshouldBeThereFN);
				checkMissing(tsvThatshouldBeThere);
			}
	
			String tsvFN = writeFolder + "tsvFileNames.txt";
			searchFilesInDirectores(tsvFN);
			ArrayList<String> outputFiles = FileUtils.readArrayList(tsvFN);
			
			combineFiles(outputFiles, combindedFN);// Combining files
			
			getMappingPercentages(mappingPercentagesFn, kallistoOutputFolder);
	
			if (transcriptsToGenesFN != null)
				sumTranscriptsToGenes(combindedFN, genesFN);
	
			if (mappingPercentagesFn != null && threshold!=0) {
				System.out.println("Removing samples where less than " + (threshold * 100) + "% of reads map");
				String writeFNtranscritps = FileUtils.replaceEnd(combindedFN, "_" + threshold + ".txt.gz");
				keepThresholdSamples(combindedFN, writeFNtranscritps);
				String writeFNgenes = FileUtils.replaceEnd(genesFN, "_" + threshold + ".txt.gz");
				keepThresholdSamples(genesFN, writeFNgenes);
	
				System.out.println("Removing bad samples (where less then" + minPercentageFeaturesExpressed
						+ "% of the genes/transcripts are expressed)");
				removeBadSamples(writeFNtranscritps);
				removeBadSamples(writeFNgenes);
	
				System.out.println("Combined expression files in folder: " + kallistoOutputFolder);
				System.out.println("Files written to: " + writeFolder);
			}
		}catch(Exception e){e.printStackTrace();}
		
	}

	private void getMappingPercentages(	String mappingPercentagesFn, String kallistoOutputFolder) throws FileNotFoundException, IOException
	{
		FileSearcher mappingPercentageFnSearcher = new FileSearcher();
		
		String mappingPercentagesFnsFn = FileUtils.makeFolderNameEndWithSlash(kallistoOutputFolder)+"errFilenames.txt";
		mappingPercentageFnSearcher.setWriteName(mappingPercentagesFnsFn);
		mappingPercentageFnSearcher.setSearchStrings(new String[]{".err"});	
		mappingPercentageFnSearcher.setFolders(kallistoOutputFolder);
		mappingPercentageFnSearcher.run();

		BufferedWriter mappingPercentageWriter = FileUtils.createWriter(mappingPercentagesFn);
		ArrayList<String> mappingFns = FileUtils.readArrayList(mappingPercentagesFnsFn);
		mappingPercentageWriter.write("Sample\tMappingPercentages\n");
		for(String mappingFn : mappingFns)
		{
			double mappingPercentage = getMappingpercentage(mappingFn);
			mappingPercentageWriter.write(new File(mappingFn).getName().replace(".err","")+"\t"+mappingPercentage + "\n");
		}
		mappingPercentageWriter.close();
	}

	private double getMappingpercentage(String mappingFn) throws FileNotFoundException, IOException
	{
		BufferedReader mappingReader = FileUtils.createReader(mappingFn);
		String line = null;
		while((line = mappingReader.readLine())!=null)
		{
			if(!line.startsWith("[quant] processed"))
				continue;
			String totalReads = line.split(" processed ")[1].split(" reads, ")[0].replace(",", "");
			String mappedReads = line.split(" reads, ")[1].split(" reads pseudoaligned")[0].replace(",", "");
			double mappingPercentage = Double.parseDouble(mappedReads)/Double.parseDouble(totalReads);
			return mappingPercentage;
		}
		return 0;
	}

	private void removeBadSamples(String filename) throws IOException {
		RemoveBadSamples.main(new String[] { "filename=" + filename, "writeFolder=" + writeFolder,
				"minPercentageExpressed=" + minPercentageFeaturesExpressed });

	}

	private void keepThresholdSamples(String combindedFN, String writeFN)
			throws FileNotFoundException, IOException {
		KeepThresholdSamples.main(new String[] { "filename=" + combindedFN,
				"mappingPercentageFN=" + mappingPercentagesFn, "writeFN=" + writeFN, "threshold=" + threshold });
	}

	private void sumTranscriptsToGenes(String combindedFN, String genesFN) {
		SumTranscriptsToGenes sumTranscriptsToGenes = new SumTranscriptsToGenes();
		sumTranscriptsToGenes.setCountFN(combindedFN);
		sumTranscriptsToGenes.setTranscriptToGeneFN(transcriptsToGenesFN);
		sumTranscriptsToGenes.setWriteName(genesFN);
		sumTranscriptsToGenes.run();
	}

	private void checkMissing(ArrayList<String> tsvFiles) throws IOException {
		File file = new File(kallistoOutputFolder);
		String parentPath = file.getAbsolutePath();
		String writeFileName = parentPath + "/Counts.txt";
		if (kallistoColumn == 3)
			writeFileName = writeFileName.replace(".txt", "_TPM_.txt");

		BufferedWriter missingWriter = FileUtils.createWriter(writeFolder + "missing.txt");
		for (int f = 0; f < tsvFiles.size(); f++) {
			// System.out.println("f=" + f + " " +tsvFiles.get(f));
			String line = tsvFiles.get(f);
			String[] columns = line.split("\t");
			String fileName = columns[0]+".gz";
			File fl = new File(fileName);
			if (!fl.exists())
				missingWriter.write("This file is missing:\t" + line + "\n");
		}
		missingWriter.close();
	}

	private void searchFilesInDirectores(String tsvFN) throws Exception {
		FileSearcher
				.main(new String[] { "folderName=" + kallistoOutputFolder, "searchStrings=.tsv", "writeFN=" + tsvFN });
	}

	private void combineFiles(ArrayList<String> tsvFiles, String combindedFN) // tsvFiles are the Kallisto count Output files
	{
		File file = new File(kallistoOutputFolder);
		String parentPath = file.getAbsolutePath();
		String writeFileName = parentPath + "/Counts.txt";
		if (kallistoColumn == 3)
			writeFileName = writeFileName.replace(".txt", "_TPM_.txt");

		MyMatrix output = null;
		for (int f = 0; f < tsvFiles.size(); f++) {
			if (f % 100 == 0)
				System.out.println("f=" + f + "/" + tsvFiles.size() + " = " + tsvFiles.get(f));
			String fileName = tsvFiles.get(f);
			printStatus(f, tsvFiles, fileName);

			MyMatrix counts = new MyMatrix(fileName);
			output = createMatrixIfNull(output, counts, tsvFiles);

			output = addValuesToMatrix(fileName, output, f, counts);
		}

		output.write(combindedFN);
		System.out.println("File written to: " + combindedFN);
	}

	private MyMatrix addValuesToMatrix(String fileName, MyMatrix output, int outC, MyMatrix counts) {
		File tempName = new File(fileName);
		File folder = new File(tempName.getParent());
		output.colNames[outC] = folder.getName();
		if (counts.rows() != output.rows())
			System.out.println(
					"WARNING! THIS FILE CONTAINS A DIFFERENT NUMBER OF ROWS INDICATING KALLISTO FAILED TO COMPLETE SUCCESSFULLY on:\n"
							+ counts.colNames[outC] + "\n" + "or:\n" + counts.colNames[0]);
		Hashtable<String, Integer> outputRowIndex = output.getRowHash();
		for (int x = 0; x < counts.rowNames.length; x++) {
			if (outputRowIndex.get(counts.getRowHeaders()[x]) == null)
				continue;
			// System.out.println(x + " " + counts.getRowHeaders()[x]);
			int row = outputRowIndex.get(counts.getRowHeaders()[x]);
			output.values[row][outC] = counts.values[x][kallistoColumn];
		}
		System.out.println("Samples added to matrix= " + outC);
		return output;
	}

	private static MyMatrix createMatrixIfNull(MyMatrix output, MyMatrix counts, ArrayList<String> tsvFiles) {
		if (output == null) {
			output = new MyMatrix(counts.rowNames.length, tsvFiles.size());
			output.rowNames = counts.rowNames;
		}
		return output;
	}

	private static void printStatus(int f, ArrayList<String> tsvFiles, String fileName) {
		if (f % 100 == 0) {
			System.out.println("File: " + f + "/" + tsvFiles.size());
			System.out.println(fileName);
			PCA.log(Double.toString(((double) f) / ((double) tsvFiles.size())));
		}
	}

	public String getKallistoOutputFolder()
	{
		return kallistoOutputFolder;
	}

	public void setKallistoOutputFolder(String kallistoOutputFolder)
	{
		this.kallistoOutputFolder = kallistoOutputFolder;
	}

	public String getWriteFolder()
	{
		return writeFolder;
	}

	public void setWriteFolder(String writeFolder)
	{
		this.writeFolder = writeFolder;
	}

	public int getKallistoColumn()
	{
		return kallistoColumn;
	}

	public void setKallistoColumn(int kallistoColumn)
	{
		this.kallistoColumn = kallistoColumn;
	}

	public String getTranscriptsToGenesFN()
	{
		return transcriptsToGenesFN;
	}

	public void setTranscriptsToGenesFN(String transcriptsToGenesFN)
	{
		this.transcriptsToGenesFN = transcriptsToGenesFN;
	}

	public double getThreshold()
	{
		return threshold;
	}

	public void setThreshold(double threshold)
	{
		this.threshold = threshold;
	}

	public int getMinPercentageFeaturesExpressed()
	{
		return minPercentageFeaturesExpressed;
	}

	public void setMinPercentageFeaturesExpressed(int minPercentageFeaturesExpressed)
	{
		this.minPercentageFeaturesExpressed = minPercentageFeaturesExpressed;
	}

	public String getMappingPercentagesFN()
	{
		return mappingPercentagesFn;
	}

	public void setMappingPercentagesFN(String mappingPercentagesFN)
	{
		this.mappingPercentagesFn = mappingPercentagesFN;
	}

	public String getTsvThatshouldBeThereFN()
	{
		return tsvThatshouldBeThereFN;
	}

	public void setTsvThatshouldBeThereFN(String tsvThatshouldBeThereFN)
	{
		this.tsvThatshouldBeThereFN = tsvThatshouldBeThereFN;
	}

	public String getGenesFN()
	{
		return genesFN;
	}

	public void setGenesFN(String genesFN)
	{
		this.genesFN = genesFN;
	}

//	static void checkArgs(String[] args) 
//	{
//		if (System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko") && args.length < 1)
//			return;
//		if (args.length < 1) {
//			System.out.println("Script requires the following argumetns:\n"
//					+ "1. kallistooutputfolder=<kallistoOutputRootFolder> - Folder containing all the kallisto output files (subfolders are allowed)\n"
//					+ "2. writeFolder=<writeFolderFN.txt> - Folder where the files will be written (default=parentFolder(<folder>))\n"
//					+ "3. kallistocColumn=<2> - Kallisto column to use (2 for counts, 3 for tpm)(default=2)\n"
//					+ "4. transcriptsToGenesFN=<transcriptstogenesfn.txt> - File that has transcripts names in the first column and corresponding genenames in the second\n"
//					+ "5. threshold=<0.7> - The minimum percentage/100 of reads that has to map for the sample to be included (default=0.7)\n"
//					+ "6. minpercentagefeaturesexpressed=<10> - Minimum percentage of transcrips that needs to be expressed for the sample to be included  (default=10)\n"
//					+ "7. mappingpercentagesfn=<mappingpercentagesfn.txt> - File that has the in the filenames(incl. directory) if the .err 1s column and the mapping percentages in the 2nd column\n"
//					+ "8. tsvthatshouldbetherefn=<tsvthatshouldbetherefn.txt> - File that contains the filenames(incl. direcory) of all .tsv files that should be there in the first column\n"
//					);
//			System.exit(1);
//		}
//
//		for (int a = 0; a < args.length; a++) {
//			System.out.println("args:\t" + args[a]);
//			String arg = args[a].split("=")[0];
//			String value = args[a].split("=")[1];
//			switch (arg.toLowerCase()) {
//			// var = new JSONutil<Vars>().read(var.JSON_FN, var);
//			case "kallistooutputfolder":
//				kallistoOutputFolder = value;
//				break;
//			case "writefolder":
//				writeFolder = value;
//				break;
//			case "kallistocolumn":
//				kallistoColumn = Integer.parseInt(value);
//				break;
//			case "transcriptstogenesfn":
//				transcriptsToGenesFN = value;
//				break;
//			case "threshold":
//				threshold = Double.parseDouble(value);
//				break;
//			case "minpercentagefeaturesexpressed":
//				minPercentageFeaturesExpressed = Integer.parseInt(value);
//				break;
//			case "mappingpercentagesfn":
//				mappingPercentagesFN = value;
//				break;
//			case "tsvthatshouldbetherefn":
//				tsvThatshouldBeThereFN = value;
//				break;
//			default:
//				System.out.println("Incorrect argument supplied:\n" + args[a] + "\nexiting");
//				System.exit(1);
//			}
//		}
//	}
}