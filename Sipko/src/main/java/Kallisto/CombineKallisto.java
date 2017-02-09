package Kallisto;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import JuhaPCA.PCA;
import Kallisto.FastQtoExpression.Vars;
import PCA.Matrix;
import Tools.FileUtils;
import Tools.JSONutil;
import Tools.FileSearcher;

public class CombineKallisto {

	static String kallistoOutputFolder = "E:/Groningen/Data/Juha/Genes31995/RemoveDuplicateTranscripts/Samples/IncDups/";
	static String writeFolder = null;
	static int kallistoColumn = 2;// 2 is read counts, 3 is TPM values in
									// kallisto output
	static String transcriptsToGenesFN = "E:/Groningen/Data/Juha/Genes31995/RemoveDuplicateTranscripts/DuplicateTranscriptsAndEnsemblMapping/ENSt2ENSg83.txt";// "E:/Groningen/Data/Annotation/hg19.v75.cdna.all.enst2ensg.txt";
	static double threshold = 0.7;
	static int minPercentageFeaturesExpressed = 10;// features are either
													// transcripts or genes
	static String mappingPercentagesFN = null;// "E:/Groningen/Test/JSON/ServerTest/Kallisto/mappingPerSample.txt";
	static String tsvThatshouldBeThereFN = null;// "E:/Groningen/Test/JSON/ServerTest/Kallisto/scriptNumberToFiles.txt";//should
												// have the .tsv filenames in
												// the first column (including
												// path)

	public static void main(String[] args) throws Exception {
		checkArgs(args);
		if (writeFolder == null)
			writeFolder = new File(kallistoOutputFolder).getAbsolutePath() + "/";

		if (tsvThatshouldBeThereFN != null) {
			ArrayList<String> tsvThatshouldBeThere = FileUtils.readArrayList(tsvThatshouldBeThereFN);
			checkMissing(tsvThatshouldBeThere);
		}

		String tsvFN = writeFolder + "tsvFileNames.txt";
		searchFilesInDirectores(tsvFN);
		ArrayList<String> outputFiles = FileUtils.readArrayList(tsvFN);
		String combindedFN = writeFolder + "counts.txt.gz";
		combineFiles(outputFiles, combindedFN);// Combining files
		String genesFN = FileUtils.replaceEnd(combindedFN, "_GENES.txt.gz");

		if (transcriptsToGenesFN != null)
			sumTranscriptsToGenes(combindedFN, genesFN);

		if (mappingPercentagesFN != null) {
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
	}

	private static void removeBadSamples(String filename) throws IOException {
		RemoveBadSamples.main(new String[] { "filename=" + filename, "writeFolder=" + writeFolder,
				"minPercentageExpressed=" + minPercentageFeaturesExpressed });

	}

	private static void keepThresholdSamples(String combindedFN, String writeFN)
			throws FileNotFoundException, IOException {
		KeepThresholdSamples.main(new String[] { "filename=" + combindedFN,
				"mappingPercentageFN=" + mappingPercentagesFN, "writeFN=" + writeFN, "threshold=" + threshold });
	}

	private static void sumTranscriptsToGenes(String combindedFN, String genesFN) {
		SumTranscriptsToGenes.main(new String[] { "filename=" + combindedFN,
				"transcriptToGeneFN=" + transcriptsToGenesFN, "writeFN=" + genesFN });
	}

	private static void checkMissing(ArrayList<String> tsvFiles) throws IOException {
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
			String fileName = columns[0];
			File fl = new File(fileName);
			if (!fl.exists())
				missingWriter.write("This file is missing:\t" + line + "\n");
		}
		missingWriter.close();
	}

	private static void searchFilesInDirectores(String tsvFN) throws Exception {
		FileSearcher
				.main(new String[] { "folderName=" + kallistoOutputFolder, "searchStrings=.tsv", "writeFN=" + tsvFN });
	}

	private static void combineFiles(ArrayList<String> tsvFiles, String combindedFN) // tsvFiles
																						// are
																						// the
																						// Kallisto
																						// count
																						// Output
																						// files
	{
		File file = new File(kallistoOutputFolder);
		String parentPath = file.getAbsolutePath();
		String writeFileName = parentPath + "/Counts.txt";
		if (kallistoColumn == 3)
			writeFileName = writeFileName.replace(".txt", "_TPM_.txt");

		Matrix output = null;
		for (int f = 0; f < tsvFiles.size(); f++) {
			if (f % 100 == 0)
				System.out.println("f=" + f + "/" + tsvFiles.size() + " = " + tsvFiles.get(f));
			String fileName = tsvFiles.get(f);
			printStatus(f, tsvFiles, fileName);

			Matrix counts = new Matrix(fileName);
			output = createMatrixIfNull(output, counts, tsvFiles);

			output = addValuesToMatrix(fileName, output, f, counts);
		}

		output.write(combindedFN);
		System.out.println("File written to: " + combindedFN);
	}

	private static Matrix addValuesToMatrix(String fileName, Matrix output, int outC, Matrix counts) {
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

	private static Matrix createMatrixIfNull(Matrix output, Matrix counts, ArrayList<String> tsvFiles) {
		if (output == null) {
			output = new Matrix(counts.rowNames.length, tsvFiles.size());
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

	static void checkArgs(String[] args) {
		if (System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko") && args.length < 1)
			return;
		if (args.length < 1) {
			System.out.println("Script requires the following argumetns:\n"
					+ "1. kallistooutputfolder=<kallistoOutputRootFolder> - Folder containing all the kallisto output files (subfolders are allowed)\n"
					+ "2. writeFolder=<writeFolderFN.txt> - Folder where the files will be written (default=parentFolder(<folder>))\n"
					+ "3. kallistocClumn=<2> - Kallisto column to use (2 for counts, 3 for tpm)(default=2)\n"
					+ "4. transcriptsToGenesFN=<transcriptstogenesfn.txt> - File that has transcripts names in the first column and corresponding genenames in the second\n"
					+ "5. threshold=<0.7> - The minimum percentage/100 of reads that has to map for the sample to be included (default=0.7)\n"
					+ "6. minpercentagefeaturesexpressed=<10> - Minimum percentage of transcrips that needs to be expressed for the sample to be included  (default=10)\n"
					+ "7. mappingpercentagesfn=<mappingpercentagesfn.txt> - File that has the in the filenames(incl. directory) if the .err 1s column and the mapping percentages in the 2nd column\n"
					+ "8. tsvthatshouldbetherefn=<tsvthatshouldbetherefn.txt> - File that contains the filenames(incl. direcory) of all .tsv files that should be there in the first column\n");
			System.exit(1);
		}

		for (int a = 0; a < args.length; a++) {
			System.out.println("args:\t" + args[a]);
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase()) {
			// var = new JSONutil<Vars>().read(var.JSON_FN, var);
			case "kallistooutputfolder":
				kallistoOutputFolder = value;
				break;
			case "writefolder":
				writeFolder = value;
				break;
			case "kallistocolumn":
				kallistoColumn = Integer.parseInt(value);
				break;
			case "transcriptstogenesfn":
				transcriptsToGenesFN = value;
				break;
			case "threshold":
				threshold = Double.parseDouble(value);
				break;
			case "minpercentagefeaturesexpressed":
				minPercentageFeaturesExpressed = Integer.parseInt(value);
				break;
			case "mappingpercentagesfn":
				mappingPercentagesFN = value;
				break;
			case "tsvthatshouldbetherefn":
				tsvThatshouldBeThereFN = value;
				break;
			default:
				System.out.println("Incorrect argument supplied:\n" + args[a] + "\nexiting");
				System.exit(1);
			}
		}
	}
}