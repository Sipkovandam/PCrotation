package Kallisto;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;

import STAR.STAR_ClusterHandler;
import Slurm.ClusterHandler;
import Slurm.SlurmJob;
import Tools.FileUtils;

public class Kallisto_ClusterHandler extends SlurmJob implements Cloneable, Serializable
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	private String kallistoIndexFileComment = "/root/directory/hg19.v75.cdna.all.42.2.idx; MANDATORY //Annotation file Kallisto should use. This file has to be created by Kallisto prior to running this script";
	private String kallistoIndexFile = "/groups/umcg-wijmenga/tmp04/umcg-svandam/Data/RNAseq/Annotation/hg19.v75.cdna.all.42.2.idx";
	private String extraArguments = "--fusion";

	private String pairStringsComment = "[\"_R1_\",\"_R2_\",\"_1.fq\",\"_2.fq\"]; MANDATORY for paired end data // comma separated list of string pairs defining the difference between forward read and backward read files.  For example, immagine a file name fastqFile_1.fq - to obtain the complementary file in this file name _R1_ is replaced iwth _R2_ and _1.fq is replaced with _2.fq. Since _R1_ is not present in the file name but _1.fq is, the complementary file becomes fastqFile_2.fq and these 2 files then are used by STAR. If you files do not contain any of these strings STAR will map the data as if it was single end data";
	private String[] pairStrings = new String[] { "_R1", "_R2", "_1.fq", "_2.fq" };

	@Override
	public void createSlurmFiles(	String fileNames,
									ClusterHandler slurmVars) throws Exception
	{
		ArrayList<String> fastQ_Files = FileUtils.readArrayList(fileNames);

		System.out.println("Creating slurm files in:\n" + slurmVars.getScriptsFolderName());
		int scriptNumber = 0;
		int fileNumber = 0;
		String tsvFilesToShScriptFN = new File(slurmVars.getScriptsFolderName()).getParent() + "/scriptNumberToFiles.txt";
		BufferedWriter tsvFilenameWriter = FileUtils.createWriter(tsvFilesToShScriptFN);
		BufferedWriter writer = null;
		out: for (int f = 0; f < fastQ_Files.size(); f++)
		{
			if (f == 0 || fileNumber >= slurmVars.getBatchSize())
			{
				scriptNumber++;
				String shellFN = slurmVars.getScriptsFolderName() + scriptNumber + ".sh";
				if (writer != null)
					writer.close();
				;
				writer = FileUtils.createWriter(shellFN);
				writeSlurmCommands(	writer,
									scriptNumber,
									slurmVars);

				fileNumber = 0;
			}
			String fastqFN = fastQ_Files.get(f);

			// continue if it is the second file of a paired end sequenced sample
			String pairedStringForward = null;
			String pairedStringReverse = null;
			for (int p = 1; p < getPairStrings().length; p += 2)
				if (fastqFN.contains(getPairStrings()[p]))
					continue out;

			for (int p = 0; p < getPairStrings().length; p += 2)
				if (fastqFN.contains(getPairStrings()[p]))
				{
					pairedStringForward = getPairStrings()[p];
					pairedStringReverse = getPairStrings()[p + 1];
				}

			String writeFolder = writeCommandsMapper(	writer,
														fastqFN,
														pairedStringForward,
														pairedStringReverse,
														scriptNumber,
														tsvFilenameWriter,
														slurmVars);

			writer.write("for file in `ls " + writeFolder + " | grep -v .bam`; do gzip $file; done;\n");
			fileNumber++;

		}
		if (writer != null)
			writer.close();
		tsvFilenameWriter.close();
	}

	private void writeSlurmCommands(BufferedWriter writer,
									int scriptNumber,
									ClusterHandler slurmVars) throws Exception
	{
		slurmVars.writeSlurmCommands(	writer,
										scriptNumber);
		writer.write("kallisto version\n");
	}

	private String writeCommandsMapper(	BufferedWriter writer,
										String fn,
										String pairedStringForward,
										String pairedStringReverse,
										int fileNumber,
										BufferedWriter tsvFilenameWriter,
										ClusterHandler slurmVars) throws Exception
	{

		File file = new File(fn);
		String writeFolderName = null;

		if (pairedStringForward != null && file.getName().contains(pairedStringForward))
			writeFolderName = writePairedEndCommand(file,
													pairedStringForward,
													pairedStringReverse,
													writer,
													fileNumber,
													tsvFilenameWriter,
													slurmVars);
		else if (file.getName().contains(".fq")) // single end
			writeFolderName = writeSingleEndCommand(file,
													writer,
													fileNumber,
													tsvFilenameWriter,
													slurmVars);
		return writeFolderName;
	}

	private String writePairedEndCommand(	File file,
											String pairedStringForward,
											String pairedStringReverse,
											BufferedWriter writer,
											int fileNumber,
											BufferedWriter tsvFilenameWriter,
											ClusterHandler slurmVars) throws IOException
	{
		String outputFolder = file.getName().split(pairedStringForward)[0] + "/";

		// if(iris)//Iris samples need a special treatment because the name of
		// the file is 1 folder higher in the folder hierarchie ><
		// {
		// String[] parentFolders = file.getParent().split("\\\\");
		// outputFolder = file.getName().replace("_1.fq.gz", "")+"/";
		// }
		writer.write("mkdir " + slurmVars.getSTAR_Folder() + outputFolder + "\n");
		String fileName = "\"" + file.getName() + "\"";
		fileName = fileName.replace("\\",
									"/");
		String fastqWithPath = file.getPath().replace(	"\\",
														"/");

		writeKallistoLinesPaired(	fastqWithPath,
									pairedStringForward,
									pairedStringReverse,
									outputFolder,
									writer,
									tsvFilenameWriter,
									fileNumber,
									slurmVars);

		return slurmVars.getSTAR_Folder() + outputFolder;
	}

	private void writeKallistoLinesPaired(	String fastqWithPath,
											String pairedStringForward,
											String pairedStringReverse,
											String outputFolder,
											BufferedWriter writer,
											BufferedWriter tsvFilenameWriter,
											int fileNumber,
											ClusterHandler slurmVars) throws IOException
	{
		String line = "kallisto quant --bias ";
		if(!slurmVars.isSharkCluster())
			line += "-t " + slurmVars.getThreads();
		line += "-i " + getKallistoIndexFile() + " -o " + slurmVars.getSTAR_Folder() + outputFolder + " " + getExtraArguments() + " \"" + fastqWithPath + "\" \"" + fastqWithPath.replace(pairedStringForward,
																																																														pairedStringReverse)
				+ "\" &> " + slurmVars.getSTAR_Folder() + outputFolder + outputFolder.replace(	"/",
																								"")
				+ ".err" + "\n";
		writer.write(line);
		//		if(line.contains("--bias"))
		//			writer.write("gzip " + sTAR_Folder+outputFolder+"fusion.txt\n");
		tsvFilenameWriter.write(slurmVars.getSTAR_Folder() + outputFolder + "abundance.tsv" + "\t" + fileNumber + ".sh" + "\t" + fastqWithPath + "\n");
	}

	private void writeCommandsFeatureCounts(BufferedWriter writer,
											String outputFolder,
											STAR_ClusterHandler sTARv) throws IOException
	{
		String alingedName = outputFolder + "Aligned.out." + sTARv.getOutMode().toLowerCase().replace(	" unsorted",
																										"").replace(" sortedbycoordinate",
																													"");
		String featureCountsLine = sTARv.getFeatureCounts() + " " + sTARv.getFeatureCountsOptions() + " -a " + sTARv.getGTFfile() + " -o " + outputFolder + "featureCounts.out " + alingedName;
		writer.write(featureCountsLine + "\n");
		if (checkKeep(	alingedName,
						sTARv))
			writer.write("rm " + alingedName + "\n");
	}

	private boolean checkKeep(	String alingedName,
								STAR_ClusterHandler sTARv)
	{
		if (sTARv.getKeepBAMsContaining() == null)
			return true;
		for (String keepString : sTARv.getKeepBAMsContaining())
		{
			if (keepString == null)
				continue;
			if (alingedName != null && alingedName.contains(keepString))
				return false;

		}
		return true;
	}

	private String writeSingleEndCommand(	File file,
											BufferedWriter writer,
											int fileNumber,
											BufferedWriter tsvFilenameWriter,
											ClusterHandler slurmvars) throws IOException
	{
		String outputFolder = file.getName().replace(	".fq.gz",
														"/");

		// if(iris)//Iris samples need a special treatment because the name of
		// the file is 1 folder higher in the folder hierarchie ><
		// {
		// String[] parentFolders = file.getParent().split("\\\\");
		// outputFolder = parentFolders[parentFolders.length-1]+"/";
		// }
		writer.write("mkdir " + slurmvars.getSTAR_Folder() + outputFolder + "\n");
		String fileName = "\"" + file.getName() + "\"";
		fileName = fileName.replace("\\",
									"/");
		String fastqWithPath = file.getPath().replace(	"\\",
														"/");

		String line = "kallisto quant --bias -t " + slurmvars.getThreads() + " --single --fragment-length=200 --sd=20 -i " + getKallistoIndexFile() + " -o " + slurmvars.getSTAR_Folder() + outputFolder + " " + getExtraArguments() + " \"" + fastqWithPath + "\" &> " + slurmvars.getSTAR_Folder() + outputFolder + "kallisto_" + outputFolder.replace(	"/",
																																																																																							"")
				+ ".err" + "\n";
		writer.write(line);
		tsvFilenameWriter.write(slurmvars.getSTAR_Folder() + outputFolder + "abundance.tsv" + "\t" + fileNumber + ".sh" + "\t" + fastqWithPath + "\n");

		return outputFolder;
	}

	public String getKallistoIndexFile()
	{
		return kallistoIndexFile;
	}

	public void setKallistoIndexFile(String kallistoIndexFile)
	{
		this.kallistoIndexFile = kallistoIndexFile;
	}

	public String getExtraArguments()
	{
		return extraArguments;
	}

	public void setExtraArguments(String extraArguments)
	{
		this.extraArguments = extraArguments;
	}

	public String[] getPairStrings()
	{
		return pairStrings;
	}

	public void setPairStrings(String[] pairStrings)
	{
		this.pairStrings = pairStrings;
	}

}
