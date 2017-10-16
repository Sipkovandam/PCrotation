package STAR;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import Slurm.ClusterHandler;
import Slurm.SlurmJob;
import Tools.FileUtils;

public class FastqDeduplicator_ClusterHandler extends SlurmJob
{
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
			String fastqFn = fastQ_Files.get(f);

			// skip to next, if it is the second file of a paired end sequenced sample
			String pairedStringForwardFastqFn = fastqFn;
			String pairedStringReverseFastqFn = null;
			for (int p = 1; p < getPairStrings().length; p += 2)
				if (fastqFn.contains(getPairStrings()[p]))
					continue out;

			for (int p = 0; p < getPairStrings().length; p += 2)
				if (fastqFn.contains(getPairStrings()[p]))
				{
					pairedStringForwardFastqFn = fastqFn;
					pairedStringReverseFastqFn = fastqFn.replace(	getPairStrings()[p],
					                                              	getPairStrings()[p + 1]);
				}

			String writeFolder = writeCommandsMapper(	writer,
			                                         	pairedStringForwardFastqFn,
														pairedStringReverseFastqFn,
														scriptNumber,
														tsvFilenameWriter,
														slurmVars);

			fileNumber++;

		}
		if (writer != null)
			writer.close();
		;
		tsvFilenameWriter.close();

	}

	private void writeSlurmCommands(BufferedWriter writer,
									int scriptNumber,
									ClusterHandler slurmVars) throws Exception
	{
		slurmVars.writeSlurmCommands(	writer,
										scriptNumber);
	}

	public String writeCommandsMapper(	BufferedWriter writer,
										String pairedStringForward,
										String pairedStringReverse,
										int fileNumber,
										BufferedWriter tsvFilenameWriter,
										ClusterHandler slurmVars) throws Exception
	{
		String writeFolderName = null;

		if (pairedStringReverse != null)
			writeFolderName = writeBBmapLines_PairedEnd(pairedStringForward,
														pairedStringReverse,
														writer,
														tsvFilenameWriter,
														slurmVars);

		else if (pairedStringForward.contains(".fq") || pairedStringForward.contains(".fastq")) // single end
			writeFolderName = writeBBmapLines_SingleEnd(pairedStringForward,
														writer,
														tsvFilenameWriter,
														slurmVars);
		return writeFolderName;
	}

	private String writeBBmapLines_SingleEnd(	String fastqFn,
											BufferedWriter writer,
											BufferedWriter tsvFilenameWriter,
											ClusterHandler slurmvars) throws IOException
	{
		String line = "dedupe.sh usejni=t";
		line += " t=" + slurmvars.getThreads();
		line += " in=" + fastqFn;
		
		String tempName = FileUtils.removeExtention(fastqFn) + "_temp.fq";
		line += " out=" + tempName;
		writer.write(line + "\n");
		
		writer.write("rm " + tempName + "\n");
		writer.write("rm " + fastqFn + "\n");
		return new File(fastqFn).getParent();
	}

	private String writeBBmapLines_PairedEnd(	String pairedStringForward,
												String pairedStringReverse,
												BufferedWriter writer,
												BufferedWriter tsvFilenameWriter,
												ClusterHandler slurmvars) throws IOException
	{
		String line = "dedupe.sh usejni=t";

		line += " t=" + slurmvars.getThreads();
		line += " in1=" + pairedStringForward;
		line += " in2=" + pairedStringReverse;

		String tempName = FileUtils.removeExtention(pairedStringForward) + "_temp.fq";
		line += " out=" + tempName;
		writer.write(line + "\n");

		String reformatLine = "reformat.sh " + tempName;
		reformatLine += " out1=" + FileUtils.addBeforeExtention(	pairedStringForward,
																	"_Deduplicated");
		reformatLine += " out2=" + FileUtils.addBeforeExtention(	pairedStringReverse,
																	"_Deduplicated");
		writer.write(reformatLine + "\n");
		
		writer.write("rm " + tempName + "\n");
		writer.write("rm " + pairedStringForward + "\n");
		writer.write("rm " + pairedStringReverse + "\n");
		return new File(pairedStringForward).getParent();
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
