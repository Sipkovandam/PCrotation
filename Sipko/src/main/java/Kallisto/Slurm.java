package Kallisto;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.stream.Stream;

import org.junit.Rule;
import org.junit.rules.TemporaryFolder;

import Kallisto.FastQtoExpression.Vars;
import Tools.ExecCommand;
import Tools.FileUtils;
import Tools.JSONutil;
import Tools.Util;

public class Slurm {

	@Rule
	public TemporaryFolder tempFolder = new TemporaryFolder();

	String fastQFiles = "E:/Groningen/Test/JSON/fastQfiles.txt";
	String writeFolder = null;
	String scriptsFolderName = null;
	String logsFolder = null;
	String errorsFolder = null;
	String mapper = "Kallisto/0.42.2.1-goolf-1.7.20";// "STAR/2.5.1b-foss-2015b";
	int batchSize = 1;
	// kallisto specific arguments
	String kallistoIndexFile = "/groups/umcg-wijmenga/tmp04/umcg-svandam/Data/RNAseq/Annotation/hg19.v75.cdna.all.42.2.idx";
	// STAR specific arguments
	String genomeDir = "";
	String buildDir = null;
	String maxMemoryGenomeBuild = "180gb";
	String gTFfile = null;
	String sjdbFileChrStartEnd = null;
	String outMode = "BAM Unsorted";
	String[] keepBAMsContaining = new String[]{"RNA14-00231_S1","RNA14-00254_S7","RNA14-00258_S4"};//if a bam file contains any of these strings it is kept
	String STARguments = null;
	boolean countExpression= true;
	//featureCounts arguments
	private String featureCounts = null;
	private String featureCountsOptions = null;

	int threads = 4;
	String walltime = "05:59:00";
	String maxMemory = "8gb";
	String outputRoot = null;
	String mappingPercentagesFN = null;

	String slurmUserName = "umcg-svandam";
	String finishedemailaddress = "sipkovandam@gmail.com";
	String tsvFilesToShScriptFN = null;

	// pairStrings to be in pairs. String in uneven index is replaced with that
	// in the even index to get the paired end partner.
	// All file names containing a string in one of the uneven indexes is
	// skipped)
	String[] pairStrings = new String[] { "_R1_", "_R2_", "_1.fq", "_2.fq" };

	public void run(String[] args) throws Exception {
		checkArgs(args);
		if (writeFolder == null)
			writeFolder = new File(fastQFiles).getParent() + "/";
		System.out.println("writeFolder=" +writeFolder);
		System.out.println("scriptsFolderName=" + scriptsFolderName);
		scriptsFolderName = writeFolder + "Scripts/";
		logsFolder = writeFolder + "Logs/";
		errorsFolder = writeFolder + "Errors/";
		outputRoot = writeFolder + "Results/";
		mappingPercentagesFN = new File(outputRoot).getParent() + "/mappingPerSample.txt";
		if (tsvFilesToShScriptFN == null)
			tsvFilesToShScriptFN = new File(scriptsFolderName).getParent() + "/scriptNumberToFiles.txt";
		System.out.println("scriptsFolderName2=" + scriptsFolderName);

		ArrayList<String> fastQ_Files = FileUtils.readArrayList(fastQFiles);
		FileUtils.makeDir(scriptsFolderName);
		FileUtils.makeDir(logsFolder);
		FileUtils.makeDir(errorsFolder);
		FileUtils.makeDir(outputRoot);
		System.out.println("scriptsFolderName=" + scriptsFolderName);
		System.out.println("countExpression=" + countExpression);
		createSlurmFiles(fastQ_Files);

		runSlurmScripts();
		sjdbFileChrStartEnd=null;
	}

	private void runSlurmScripts() throws FileNotFoundException, IOException {
		// carefull with this. Kallisto200sh script does not update! You need to
		// delete and reimport project of that!
		tempFolder.create();
		
		String kallistoShellFN = "Kallisto200.sh";//FileUtils.prepareBinaryFromJar("Kallisto200.sh");
		
		System.out.println("exists = " + new File(kallistoShellFN).exists()+ "\t" + kallistoShellFN);
		String command = "bash " + kallistoShellFN + " " + scriptsFolderName + "*.sh " + outputRoot + " "
				+ slurmUserName + " " + mappingPercentagesFN + " " + finishedemailaddress;
		System.out.println("Shellcommand = " + command);
		ExecCommand exec = new ExecCommand(command);
		// System.out.println("execute output: \n"+exec.getOutput());
		if(exec.getError().length()>1)
			System.out.println("execute error: \n" + exec.getError());
		else
			System.out.println("No errors occurred running the slurm scripts");
	}
	// Redundant, Jesse heeft hier voor een paar classes gegeven
	// private void runShell(String command) {
	// try
	// {
	// Runtime rt = Runtime.getRuntime();
	// Process proc = rt.exec(command);
	// InputStream stderr = proc.getErrorStream();
	// InputStreamReader isr = new InputStreamReader(stderr);
	// BufferedReader br = new BufferedReader(isr);
	// String line = null;
	// while ( (line = br.readLine()) != null)
	// System.out.println(line);
	// int exitVal = proc.waitFor();
	// System.out.println("Process exitValue: " + exitVal);
	// } catch (Throwable t)
	// {
	// t.printStackTrace();
	// }
	// }

	private void createSlurmFiles(ArrayList<String> fastQ_Files) throws Exception {
		int scriptNumber = 1;
		int fileNumber = 0;
		BufferedWriter tsvFilenameWriter = FileUtils.createWriter(tsvFilesToShScriptFN);
		BufferedWriter writer = null;
		out: for (int f = 0; f < fastQ_Files.size(); f++) {
			if(f==0 || fileNumber>=batchSize)
			{
				String shellFN = scriptsFolderName + scriptNumber + ".sh";
				if(writer !=null)
					closeWriter(writer);
				writer = FileUtils.createWriter(shellFN);
				writeSlurmCommands(writer, scriptNumber);
				scriptNumber++;
				fileNumber = 0;
			}
			String fastqFN = fastQ_Files.get(f);

			// continue if it is the second file of a paired end sequenced
			// sample
			String pairedStringForward = null;
			String pairedStringReverse = null;
			for (int p = 1; p < pairStrings.length; p += 2)
				if (fastqFN.contains(pairStrings[p]))
					continue out;

			for (int p = 0; p < pairStrings.length; p += 2)
				if (fastqFN.contains(pairStrings[p])) {
					pairedStringForward = pairStrings[p];
					pairedStringReverse = pairStrings[p + 1];
				}

			String outputFolder=writeCommandsMapper(writer, fastqFN, pairedStringForward, pairedStringReverse, scriptNumber, tsvFilenameWriter);
			if(!outMode.toLowerCase().equals("none"))
				writeCommandsFeatureCounts(writer, outputFolder);
			
			fileNumber++;
			if(buildDir!=null)
				break;
		}
		if(writer !=null)
			closeWriter(writer);
		tsvFilenameWriter.close();
	}

	private void writeCommandsFeatureCounts(BufferedWriter writer, String outputFolder) throws IOException {
		String alingedName = outputFolder+"Aligned.out."+outMode.toLowerCase().replace(" unsorted", "").replace(" sortedbycoordinate", "");
		String featureCountsLine = featureCounts+ " "+featureCountsOptions + " -a "+ gTFfile + " -o "+ outputFolder +"featureCounts.out "+alingedName;
		writer.write(featureCountsLine+"\n");
		if(checkKeep(alingedName))
			writer.write("rm "+alingedName+"\n");
	}

	private boolean checkKeep(String alingedName) {
		for(String keepString : keepBAMsContaining)
			if(alingedName.contains(keepString))
				return false;
		return true;
	}

	private void closeWriter(BufferedWriter writer) throws IOException {
		if(mapper.toLowerCase().contains("star"))
			writer.write("STAR --genomeLoad Remove --genomeDir " + genomeDir +"\n");
		writer.close();
	}

	private void writeSlurmCommands(BufferedWriter writer, int fileNumber) throws Exception {
		writer.write("#!/bin/bash" + "\n");
		writer.write("#SBATCH --job-name="+mapper + fileNumber + "\n");
		writer.write("#SBATCH --output=" + logsFolder + fileNumber + ".out\n");
		writer.write("#SBATCH --error=" + errorsFolder + fileNumber + ".err\n");
		// writer.write("#SBATCH --partition=medium\n");
		writer.write("#SBATCH --time=" + walltime + "\n");
		writer.write("#SBATCH --cpus-per-task " + threads + "\n");
		if(buildDir!=null)
			writer.write("#SBATCH --mem " + maxMemoryGenomeBuild + "\n");
		else
			writer.write("#SBATCH --mem " + maxMemory + "\n");
		writer.write("#SBATCH --nodes 1\n");
		writer.write("#SBATCH --export=NONE\n");
		writer.write("#SBATCH --get-user-env=L\n");
		// writer.write("#SBATCH --open-mode=append\n");
		
		writer.write("module load " + mapper + "\n");
		if(mapper.toLowerCase().contains("star"))
			writer.write("STAR --version\n");
		if (mapper.toLowerCase().contains("kallisto"))
			writer.write("kallisto version\n");
	}

	private String writeCommandsMapper(BufferedWriter writer, String fn, String pairedStringForward,
			String pairedStringReverse, int fileNumber, BufferedWriter tsvFilenameWriter) throws Exception {

		File file = new File(fn);
		boolean iris = false;
		String outputFolder = null;
		if (fn.contains("Iris"))// does not do anything atm (Iris changed the
								// structure of her folders)
			iris = true;
		if (pairedStringForward != null && file.getName().contains(pairedStringForward))
			outputFolder=writePairedEndCommand(file, pairedStringForward, pairedStringReverse, iris, writer, fileNumber,
					tsvFilenameWriter);
		else if (file.getName().contains(".fq")) // single end
			outputFolder=writeSingleEndCommand(file, iris, writer, fileNumber, tsvFilenameWriter);
		return outputFolder;
	}

	private String writePairedEndCommand(File file, String pairedStringForward, String pairedStringReverse,
			boolean iris, BufferedWriter writer, int fileNumber, BufferedWriter tsvFilenameWriter) throws IOException {
		String outputFolder = file.getName().split(pairedStringForward)[0] + "/";
		outputFolder=checkBuild(outputFolder);
		
		// if(iris)//Iris samples need a special treatment because the name of
		// the file is 1 folder higher in the folder hierarchie ><
		// {
		// String[] parentFolders = file.getParent().split("\\\\");
		// outputFolder = file.getName().replace("_1.fq.gz", "")+"/";
		// }
		writer.write("mkdir " + outputRoot + outputFolder + "\n");
		String fileName = "\"" + file.getName() + "\"";
		fileName = fileName.replace("\\", "/");
		String fastqWithPath = file.getPath().replace("\\", "/");
		if (mapper.toLowerCase().contains("kallisto"))
			writeKallistoLinesPaired(fastqWithPath, pairedStringForward, pairedStringReverse, outputFolder, writer,
					tsvFilenameWriter, fileNumber);
		if (mapper.toLowerCase().contains("star"))
			writeSTAR_Lines(fastqWithPath, pairedStringForward, pairedStringReverse, outputFolder, writer,
					tsvFilenameWriter, fileNumber);
		return outputRoot+outputFolder;
	}

	private void writeKallistoLinesPaired(String fastqWithPath, String pairedStringForward,
			String pairedStringReverse, String outputFolder, BufferedWriter writer, BufferedWriter tsvFilenameWriter,
			int fileNumber) throws IOException {
		String line = "kallisto quant --bias -t " + threads + " -i " + kallistoIndexFile + " -o " + outputRoot
				+ outputFolder + " \"" + fastqWithPath + "\" \""
				+ fastqWithPath.replace(pairedStringForward, pairedStringReverse) + "\" &> " + outputRoot + outputFolder
				+ outputFolder.replace("/", "") + ".err" + "\n";
		writer.write(line);
		tsvFilenameWriter.write(
				outputRoot + outputFolder + "abundance.tsv" + "\t" + fileNumber + ".sh" + "\t" + fastqWithPath + "\n");
	}

	private void writeSTAR_Lines(String fastqWithPath, 
										  String outputFolder,
										  BufferedWriter writer,
										  BufferedWriter tsvFilenameWriter,
										  int fileNumber) throws IOException 
	{
		writeSTAR_Lines(fastqWithPath, 
						null,
						null,
						outputFolder,
						writer,
						tsvFilenameWriter,
						fileNumber);
	}
	
	private void writeSTAR_Lines(String fastqWithPath, 
											  String pairedStringForward,
											  String pairedStringReverse,
											  String outputFolder,
											  BufferedWriter writer,
											  BufferedWriter tsvFilenameWriter,
											  int fileNumber) throws IOException {
		String line = "STAR ";
		if(sjdbFileChrStartEnd ==null)
			line = line + "--genomeLoad LoadAndKeep ";//The --genomeLoad LoadAndKeep option loads the genome if it's not already loaded, and then starts mapping
			line = line	+ "--genomeDir " + genomeDir +
					  " --readFilesIn " + fastqWithPath;
		if(pairedStringForward!=null)		
			line = line + " " + fastqWithPath.replace(pairedStringForward, pairedStringReverse);
		line = line	+ " --outFilterMultimapNmax 1 ";
		if(gTFfile != null && sjdbFileChrStartEnd!=null)//genome should be build with the GTF file included (See STAR manual for how)
			line = line	+ "--sjdbGTFfile " + gTFfile;
		line = line	+ " --outSAMtype " + outMode
				    + " --runThreadN " + threads 
				    + " --readFilesCommand zcat" 
				    + " --outFileNamePrefix "+ outputRoot + outputFolder;
		if(sjdbFileChrStartEnd!=null && sjdbFileChrStartEnd.length()>0 && !sjdbFileChrStartEnd.contentEquals("null"))
			line+=" --sjdbFileChrStartEnd "+sjdbFileChrStartEnd;
		if(STARguments!=null)
			line+=" "+STARguments;
		if(countExpression)
			line+=" --quantMode GeneCounts";
		writer.write(line+"\n");
		tsvFilenameWriter.write(outputRoot + outputFolder + "Aligned.out.bam" + "\t" + fileNumber + ".sh" + "\t"
				+ fastqWithPath + "\n");
	}

	private String writeSingleEndCommand(File file, boolean iris, BufferedWriter writer, int fileNumber,
			BufferedWriter tsvFilenameWriter) throws IOException {
		String outputFolder = file.getName().replace(".fq.gz", "/");
		outputFolder=checkBuild(outputFolder);
		
		// if(iris)//Iris samples need a special treatment because the name of
		// the file is 1 folder higher in the folder hierarchie ><
		// {
		// String[] parentFolders = file.getParent().split("\\\\");
		// outputFolder = parentFolders[parentFolders.length-1]+"/";
		// }
		writer.write("mkdir " + outputRoot + outputFolder + "\n");
		String fileName = "\"" + file.getName() + "\"";
		fileName = fileName.replace("\\", "/");
		String fastqWithPath = file.getPath().replace("\\", "/");

		if (mapper.toLowerCase().contains("kallisto")) {
			String line = "kallisto quant --bias -t " + threads + " --single --fragment-length=200 --sd=20 -i "
					+ kallistoIndexFile + " -o " + outputRoot + outputFolder + " \"" + fastqWithPath + "\" &> "
					+ outputRoot + outputFolder + "kallisto_" + outputFolder.replace("/", "") + ".err" + "\n";
			writer.write(line);
			tsvFilenameWriter.write(outputRoot + outputFolder + "abundance.tsv" + "\t" + fileNumber + ".sh" + "\t"
					+ fastqWithPath + "\n");
		}
		if (mapper.toLowerCase().contains("star")) {
			writeSTAR_Lines(fastqWithPath, outputFolder, writer, tsvFilenameWriter, fileNumber);
		}
		return outputFolder;
	}

	private String checkBuild(String outputFolder) {
		if(buildDir != null)
			return "";
		return outputFolder;
	}

	void checkArgs(String[] args) {
		if (System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko") && args.length < 1)
			return;
		if (args.length < 1) {
			System.out.println("Script requires the following argumetns:\n"
					+ "1.  fastQFiles=<fastQFiles.txt> - File containing all the fastQfileNames preferably including their directories\n"
					+ "2.  kallistoIndexFile=<kallistoIndexFile.idx> - Kallist index file that should be used\n"
					+ "3.  writefolder=<writefolder> - Folder where the output will be written (default=<fastQFilesRoot>/Kallisto/Scripts/)\n"
					+ "4.  kallistoVersion=<kallistoVersion> - Kallisto version that will be used (default=Kallisto/0.42.2.1-goolf-1.7.20)\n"
					+ "5.  kallistoThreads=<kallistoThreads> - Kallisto version that will be used (default=4)\n"
					+ "6.  kallistoWalltime=<kallistoWalltime> - Kallisto version that will be used (default=05:59:00)\n"
					+ "7. kallistoMaxMemory=<kallistoMaxMemory> - Kallisto version that will be used (default=8gb)\n"
					+ "8. slurmUserName=<umcg-user> - slurmUserName... (default=8gb)\n"
					+ "9. finishedemailaddress=<sipkovandam@gmail.com> - Will send an email here once finished (default=8gb)\n"
					+ "10. pairStrings=<_R1_,_R2_,_1.fq,_2.fq>) - list of string pairs that indicate differences between paired files that the script should expect\n"
					+ "PairStrings to be in pairs and are split by commas, multiple different pairs are epxected they can also be concatinated with commas behind the first pair\n"
					+ "String in uneven index is replaced with that in the even index to get the paired end partner\n"
					+ "All file names containing a string in one of the uneven indexes is skipped)");
			System.exit(1);
		}

		for (int a = 0; a < args.length; a++) {
			System.out.println("arg=" + args[a]);
			String[] split =args[a].split("=");
			String arg = split[0];
			if(split.length <2)
			{
				System.out.println("Argument missing, skipping:" + args[a]);
				continue;
			}
			String value = split[1];
			
			switch (arg.toLowerCase()) {
			// case "fastQFiles":
			// var.JSON_FN =value;
			// var = new JSONutil<Vars>().read(var.JSON_FN, var);
			// break;
			case "fastqfiles":
				fastQFiles = parseString(value);
				break;
			case "kallistoindexfile":
				kallistoIndexFile = parseString(value);
				break;
//			case "logsfolder":
//				logsFolder = parseString(value);
//				break;
//			case "errorsfolder":
//				errorsFolder = parseString(value);
//				break;
			case "kallistooutputroot":
				outputRoot = parseString(value);
				break;
			case "kallistoversion":
				mapper = parseString(value);
				break;
			case "kallistothreads":
				threads = Integer.parseInt(value);
				break;
			case "kallistowalltime":
				walltime = parseString(value);
				break;
//			case "outputroot":
//				outputRoot = parseString(value);
//				break;
			case "mapper":
				mapper = parseString(value);
				break;
			case "threads":
				threads = Integer.parseInt(value);
				break;
			case "walltime":
				walltime = parseString(value);
				break;
			case "kallistomaxmemory":
				maxMemory = parseString(value);
				break;
			case "maxmemory":
				maxMemory = parseString(value);
				break;
			case "pairstrings":
				pairStrings = value.split(",");
				break;
			case "slurmusername":
				slurmUserName = parseString(value);
				break;
			case "finishedemailaddress":
				finishedemailaddress = parseString(value);
				break;
			case "writefolder":
				writeFolder = parseString(value);
				break;
//			case "mappingpercentagesfn":
//				mappingPercentagesFN = parseString(value);
//				break;
			case "tsvfilestoshscriptfn":
				tsvFilesToShScriptFN = parseString(value);
				break;
			case "genomedir":
				genomeDir = parseString(value);
				break;
			case "gtffile":
				gTFfile = parseString(value);
				break;
			case "sjdbfilechrstartend":
				sjdbFileChrStartEnd = parseString(value);
				break;
			case "outmode":
				outMode = parseString(value);
				break;
			case "countexpression":
				countExpression = Boolean.parseBoolean(value);
				break;
			case "batchsize":
				batchSize = Integer.parseInt(value);
				break;	
			case "star_arguments":
				STARguments = parseString(value);
				break;		
			case "starguments":
				STARguments = parseString(value);
				break;
			case "builddir":
				buildDir = parseString(value);
				break;
			case "featurecounts":
				featureCounts = parseString(value);
				break;
			case "featurecountsoptions":
				featureCountsOptions = parseString(value);
				break;
			case "keepbamscontaining":
				if(parseString(value) !=null)
					keepBAMsContaining = parseString(value).split(",");
				break;
			case "maxmemorygenomebuild":
				if(parseString(value) !=null)
					maxMemoryGenomeBuild = parseString(value);
				break;
			default:
				System.out.println("Incorrect argument supplied:\n" + args[a] + "\nexiting");
				System.exit(1);
			}
		}
	}

	private String parseString(String value) {
		if(value.equals("null"))
			return null;
		return value;
	}
}
