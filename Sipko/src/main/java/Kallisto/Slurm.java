package Kallisto;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.stream.Stream;

import org.junit.Rule;
import org.junit.rules.TemporaryFolder;

import Kallisto.FastQtoExpression.Vars;
import STAR.STAR_Variables;
import Tools.ExecCommand;
import Tools.FileUtils;
import Tools.JSONutil;
import Tools.Script;
import Tools.Util;

public class Slurm <M> extends Script<Slurm<M>> implements Cloneable{

	//@Rule
	//public TemporaryFolder tempFolder = new TemporaryFolder();

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private String fastQFiles = null;//semitransient
	private String outputRoot = null;//semitransient
	private transient String scriptsFolderName = null;
	private transient String logsFolder = null;
	private transient String errorsFolder = null;
	private transient String sTAR_Folder = null;
	private String mapperComment = "STAR/2.5.1b-foss-2015b or Kallisto/0.42.2.1-goolf-1.7.20;MANDATORY // The executeable of the STAR version to be used";
	private String mapper = "STAR/2.5.1b-foss-2015b";// "STAR/2.5.1b-foss-2015b";
	private String batchSizeComment = "MANDATORY //number of files to run in 1 batch. Each batch only loads the genome once to save runtime. However if the number is to large the walltime may be hit and some files will not/partly be included in the analysis; no warnings are generated when this happens but if the walltime has been hit can be checked in the /error folder";
	private int batchSize = 1;
	//mapper variables
	public M mapperVars = null;
	//Slurm Variables
	private String threadsComment = "24; MANDATORY // Number of threads each node will use";
	private int threads = 24;
	private String walltimeComment = "23:59:00; MANDATORY // Maximum time each job has to complete the batch defined in (batchSize)";
	private String walltime = "23:59:00";
	private String maxMemoryComment = "45gb; MANDATORY // Maximum amount of memory the node can use when analysing the samples";
	private String maxMemory = "45gb";

	private String slurmUserNameComment = "umcg-svandam; MANDATORY // your slurm username";
	private String slurmUserName = "umcg-svandam";
	private String finishedemailaddressComment = "sipkovandam@gmail.com; OPTIONAL // best left empty, but sends a message whenever the mapping is completed. Both after first and second pass.";
	private String finishedemailaddress = "sipkovandam@gmail.com";

	// pairStrings to be in pairs. String in uneven index is replaced with that
	// in the even index to get the paired end partner.
	// All file names containing a string in one of the uneven indexes is
	// skipped)
	private String pairStringsComment = "[\"_R1_\",\"_R2_\",\"_1.fq\",\"_2.fq\"]; MANDATORY for paired end data // comma separated list of string pairs defining the difference between forward read and backward read files.  For example, immagine a file name fastqFile_1.fq - to obtain the complementary file in this file name _R1_ is replaced iwth _R2_ and _1.fq is replaced with _2.fq. Since _R1_ is not present in the file name but _1.fq is, the complementary file becomes fastqFile_2.fq and these 2 files then are used by STAR. If you files do not contain any of these strings STAR will map the data as if it was single end data";
	private String[] pairStrings = new String[] { "_R1_", "_R2_", "_1.fq", "_2.fq" };

	public void run() {
		try
		{
			//checkArgs(args);
			init();
	
			ArrayList<String> fastQ_Files = FileUtils.readArrayList(fastQFiles);
			createSlurmFiles(fastQ_Files);
	
			runSlurmScripts();
		}catch(Exception e){e.printStackTrace();}
	}
	
	private void init() {
		if (outputRoot == null)
			outputRoot = new File(fastQFiles).getParent() + "/";
		p("WriteFolder=" +outputRoot);
		
		scriptsFolderName = outputRoot + "Scripts/";
		logsFolder = outputRoot + "Logs/";
		errorsFolder = outputRoot + "Errors/";
		if(sTAR_Folder==null)
			sTAR_Folder = outputRoot + "Results/";
		
		if(mapperVars== null)
		{
			if(isSTAR())
				mapperVars = (M) new STAR_Variables();
			if(isKallisto())
				mapperVars = (M) new Kallisto_Variables();
			else
			{
				p("Please select an appropriate mapper containing either \"STAR\" of \"kallisto\" in the mapper name");
			}
		}

		
		FileUtils.makeDir(scriptsFolderName);
		FileUtils.makeDir(logsFolder);
		FileUtils.makeDir(errorsFolder);
		FileUtils.makeDir(sTAR_Folder);
		this.writeConfig();
	}

	private boolean isSTAR() {
		return mapper.toLowerCase().contains("star");
	}
	private boolean isKallisto() {
		return mapper.toLowerCase().contains("kallisto");
	}

	public void writeConfig(String jsonFN)
	{
		String kallistoFN = FileUtils.addBeforeExtention(jsonFN, "_kallisto");
		this.mapperVars= (M) new Kallisto_Variables();
		super.writeConfig(kallistoFN, this);
		String STAR_FN = FileUtils.addBeforeExtention(jsonFN, "_STAR");
		this.mapperVars= (M) new STAR_Variables();
		super.writeConfig(STAR_FN, this);
	}
	
	private void runSlurmScripts() throws FileNotFoundException, IOException {
		// carefull with this. Kallisto200sh script does not update! You need to
		// delete and reimport project of that!
		//tempFolder.create();
		
		String kallistoShellFN = "Kallisto200.sh";//FileUtils.prepareBinaryFromJar("Kallisto200.sh");
		
		if(!new File(kallistoShellFN).exists())
			p("Slurmscript does not exist, please copy to:\t" + kallistoShellFN);
		
		String mappingPercentagesFN = new File(sTAR_Folder).getParent() + "/mappingPerSample.txt";
		String command = "bash " + kallistoShellFN + " " + scriptsFolderName + "*.sh " + sTAR_Folder + " "
				+ slurmUserName + " " + mappingPercentagesFN + " " + finishedemailaddress;
		p("Output folder of slurm scripts:\t" + this.getOutputRoot());
		p("Running SLURM scripts using:\t" + command);
		ExecCommand exec = new ExecCommand(command);
		// p("execute output: \n"+exec.getOutput());
		if(exec.getError().length()>1)
			p("execute error: \n" + exec.getError());
		else
			p("No errors occurred running the slurm scripts");
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
	// p(line);
	// int exitVal = proc.waitFor();
	// p("Process exitValue: " + exitVal);
	// } catch (Throwable t)
	// {
	// t.printStackTrace();
	// }
	// }

	private void createSlurmFiles(ArrayList<String> fastQ_Files) throws Exception {
		p("Creating slurm files in:\n"+ scriptsFolderName);
		int scriptNumber = 0;
		int fileNumber = 0;
		String tsvFilesToShScriptFN = new File(scriptsFolderName).getParent() + "/scriptNumberToFiles.txt";
		BufferedWriter tsvFilenameWriter = FileUtils.createWriter(tsvFilesToShScriptFN);
		BufferedWriter writer = null;
		out: for (int f = 0; f < fastQ_Files.size(); f++) {
			if(f==0 || fileNumber>=batchSize)
			{
				scriptNumber++;
				String shellFN = scriptsFolderName + scriptNumber + ".sh";
				if(writer !=null)
					closeWriter(writer);
				writer = FileUtils.createWriter(shellFN);
				writeSlurmCommands(writer, scriptNumber);
				
				fileNumber = 0;
			}
			String fastqFN = fastQ_Files.get(f);

			// continue if it is the second file of a paired end sequenced sample
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

			String writeFolder = writeCommandsMapper(writer, fastqFN, pairedStringForward, pairedStringReverse, scriptNumber, tsvFilenameWriter);

			if(isGenomeBuildRun())
				break;
			writer.write("gzip " + writeFolder+"*\n");
			fileNumber++;
		}
		if(writer !=null)
			closeWriter(writer);
		tsvFilenameWriter.close();
	}

	private boolean isGenomeBuildRun() {
		if(isSTAR())
		{			
			if(this.getSTAR_V().isSaveGenome())
			return true;
		}
		return false;
	}

	private STAR_Variables getSTAR_V() {
		return (STAR_Variables) mapperVars;
	}
	public String getOutputRoot() {
		return FileUtils.makeFolderNameEndWithSlash(outputRoot);
	}

	public void setOutputRoot(String outputRoot) {
		this.outputRoot = outputRoot;
	}

	private Kallisto_Variables getKallistoV() {
		return (Kallisto_Variables) mapperVars;
	}

	private void writeCommandsFeatureCounts(BufferedWriter writer, String outputFolder, STAR_Variables sTARv) throws IOException {
		String alingedName = outputFolder+"Aligned.out."+sTARv.getOutMode().toLowerCase().replace(" unsorted", "").replace(" sortedbycoordinate", "");
		String featureCountsLine = sTARv.getFeatureCounts()+ " "+sTARv.getFeatureCountsOptions() + " -a "+ sTARv.getGTFfile() + " -o "+ outputFolder +"featureCounts.out "+alingedName;
		writer.write(featureCountsLine+"\n");
		if(checkKeep(alingedName, sTARv))
			writer.write("rm "+alingedName+"\n");
	}
	
	private boolean checkKeep(String alingedName, STAR_Variables sTARv) {
		if(sTARv.getKeepBAMsContaining()==null)
			return true;
		for(String keepString : sTARv.getKeepBAMsContaining())
			if(alingedName!= null && alingedName.contains(keepString))
				return false;
		return true;
	}

	private void closeWriter(BufferedWriter writer) throws IOException {
		if(isSTAR())
			writer.write("STAR --genomeLoad Remove --genomeDir " + getSTAR_V().getGenomeDir() +"\n");
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
		if(isGenomeBuildRun())
			writer.write("#SBATCH --mem " + getSTAR_V().getMaxMemoryGenomeBuild() + "\n");
		else
			writer.write("#SBATCH --mem " + maxMemory + "\n");
		writer.write("#SBATCH --nodes 1\n");
		writer.write("#SBATCH --export=NONE\n");
		writer.write("#SBATCH --get-user-env=L\n");
		// writer.write("#SBATCH --open-mode=append\n");
		
		writer.write("module load " + mapper + "\n");
		if(isSTAR())
			writer.write("STAR --version\n");
		if (isKallisto())
			writer.write("kallisto version\n");
	}

	private String writeCommandsMapper(BufferedWriter writer, String fn, String pairedStringForward,
			String pairedStringReverse, int fileNumber, BufferedWriter tsvFilenameWriter) throws Exception {

		File file = new File(fn);
		boolean iris = false;
		String writeFolderName = null;
		if (fn.contains("Iris"))// does not do anything atm (Iris changed the
								// structure of her folders)
			iris = true;
		if (pairedStringForward != null && file.getName().contains(pairedStringForward))
			writeFolderName=writePairedEndCommand(file, pairedStringForward, pairedStringReverse, iris, writer, fileNumber,
					tsvFilenameWriter);
		else if (file.getName().contains(".fq")) // single end
			writeFolderName=writeSingleEndCommand(file, iris, writer, fileNumber, tsvFilenameWriter);
		return writeFolderName;
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
		writer.write("mkdir " + sTAR_Folder + outputFolder + "\n");
		String fileName = "\"" + file.getName() + "\"";
		fileName = fileName.replace("\\", "/");
		String fastqWithPath = file.getPath().replace("\\", "/");
		if (isKallisto())
			writeKallistoLinesPaired(fastqWithPath, pairedStringForward, pairedStringReverse, outputFolder, writer,
					tsvFilenameWriter, fileNumber);
		if (isSTAR())
		{
			STAR_Variables STARv = getSTAR_V();
			writeSTAR_Lines(fastqWithPath, pairedStringForward, pairedStringReverse, outputFolder, writer,
					tsvFilenameWriter, fileNumber, STARv);
			
		}
		return sTAR_Folder+outputFolder;
	}

	private void writeKallistoLinesPaired(String fastqWithPath, 
										  String pairedStringForward,
										  String pairedStringReverse, 
										  String outputFolder, 
										  BufferedWriter writer, 
										  BufferedWriter tsvFilenameWriter,
										  int fileNumber) throws IOException {
		String line = "kallisto quant --bias -t " + threads + " -i " + getKallistoV().getKallistoIndexFile() + " -o " + sTAR_Folder
				+ outputFolder + " \"" + fastqWithPath + "\" \""
				+ fastqWithPath.replace(pairedStringForward, pairedStringReverse) + "\" &> " + sTAR_Folder + outputFolder
				+ outputFolder.replace("/", "") + ".err" + "\n";
		writer.write(line);
		tsvFilenameWriter.write(
				sTAR_Folder + outputFolder + "abundance.tsv" + "\t" + fileNumber + ".sh" + "\t" + fastqWithPath + "\n");
	}

	private void writeSTAR_Lines(String fastqWithPath, 
										  String outputFolder,
										  BufferedWriter writer,
										  BufferedWriter tsvFilenameWriter,
										  int fileNumber, STAR_Variables sTARv) throws IOException 
	{
		writeSTAR_Lines(fastqWithPath, 
						null,
						null,
						outputFolder,
						writer,
						tsvFilenameWriter,
						fileNumber, sTARv);
	}
	
	private void writeSTAR_Lines(String fastqWithPath, 
											  String pairedStringForward,
											  String pairedStringReverse,
											  String outputFolder,
											  BufferedWriter writer,
											  BufferedWriter tsvFilenameWriter,
											  int fileNumber, STAR_Variables sTARv) throws IOException {
		String line = "STAR ";
		if(sTARv.getSjdbFileChrStartEnd() ==null)
			line = line + "--genomeLoad LoadAndKeep ";//The --genomeLoad LoadAndKeep option loads the genome if it's not already loaded, and then starts mapping
			line = line	+ "--genomeDir " + sTARv.getGenomeDir() +
					  " --readFilesIn " + fastqWithPath;
		if(pairedStringForward!=null)		
			line = line + " " + fastqWithPath.replace(pairedStringForward, pairedStringReverse);
		if(sTARv.getGTFfile() != null && isGenomeBuildRun())//genome should be build with the GTF file included (See STAR manual for how)
			line = line	+ " --sjdbGTFfile " + sTARv.getGTFfile();
		line = line	+ " --outSAMtype " + sTARv.getOutMode()
				    + " --runThreadN " + threads 
				    + " --readFilesCommand zcat" 
				    + " --outFileNamePrefix "+ sTAR_Folder + outputFolder;
		if(sTARv.getSjdbFileChrStartEnd() !=null && sTARv.getSjdbFileChrStartEnd().length()>0 && !sTARv.getSjdbFileChrStartEnd().contentEquals("null"))
			line+=" --sjdbFileChrStartEnd "+sTARv.getSjdbFileChrStartEnd();
		if(sTARv.getSTAR_Extra_Arguments()!=null)
			line+=" "+sTARv.getSTAR_Extra_Arguments();
		if(sTARv.getCountExpression())
			line+=" --quantMode GeneCounts";
		writer.write(line+"\n");
		tsvFilenameWriter.write(sTAR_Folder + outputFolder + "Aligned.out.bam" + "\t" + fileNumber + ".sh" + "\t"
				+ fastqWithPath + "\n");
		
		if(sTARv.getOutMode() != null && !sTARv.getOutMode().toLowerCase().equals("none"))
			writeCommandsFeatureCounts(writer, outputFolder, sTARv);
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
		writer.write("mkdir " + sTAR_Folder + outputFolder + "\n");
		String fileName = "\"" + file.getName() + "\"";
		fileName = fileName.replace("\\", "/");
		String fastqWithPath = file.getPath().replace("\\", "/");

		if (isKallisto()) {
			String line = "kallisto quant --bias -t " + threads + " --single --fragment-length=200 --sd=20 -i "
					+ getKallistoV().getKallistoIndexFile() + " -o " + sTAR_Folder + outputFolder + " \"" + fastqWithPath + "\" &> "
					+ sTAR_Folder + outputFolder + "kallisto_" + outputFolder.replace("/", "") + ".err" + "\n";
			writer.write(line);
			tsvFilenameWriter.write(sTAR_Folder + outputFolder + "abundance.tsv" + "\t" + fileNumber + ".sh" + "\t"
					+ fastqWithPath + "\n");
		}
		if (isSTAR()) 
		{
			STAR_Variables STARv = getSTAR_V();
			writeSTAR_Lines(fastqWithPath, outputFolder, writer, tsvFilenameWriter, fileNumber,STARv);
		}
		return outputFolder;
	}

	private String checkBuild(String outputFolder) {
		if(isGenomeBuildRun())
			return "";
		return outputFolder;
	}

	void checkArgs(String[] args) {
		if (System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko") && args.length < 1)
			return;
		if (args.length < 1) {
			p("Script requires the following argumetns:\n"
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
			p("arg=" + args[a]);
			String[] split =args[a].split("=");
			String arg = split[0];
			if(split.length <2)
			{
				p("Argument missing, skipping:" + args[a]);
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
//			case "logsfolder":
//				logsFolder = parseString(value);
//				break;
//			case "errorsfolder":
//				errorsFolder = parseString(value);
//				break;
			case "kallistooutputroot":
				sTAR_Folder = parseString(value);
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
			case "outputroot":
				sTAR_Folder = parseString(value);
				break;
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
				outputRoot = parseString(value);
				break;
			case "batchsize":
				batchSize = Integer.parseInt(value);
				break;	
			default:
				p("Incorrect argument supplied:\n" + args[a] + "\nexiting");
				System.exit(1);
			}
		}
	}

	private String parseString(String value) {
		if(value.equals("null"))
			return null;
		return value;
	}

	public String get_STAR_Folder() {
		return sTAR_Folder;
	}

	public void set_STAR_Folder(String sTAR_Folder) {
		this.sTAR_Folder = sTAR_Folder;
	}

	public String getFastQFiles() {
		return fastQFiles;
	}

	public void setFastQFiles(String fastQFiles) {
		this.fastQFiles = fastQFiles;
	}

	public String getScriptsFolderName() {
		return scriptsFolderName;
	}

	public void setScriptsFolderName(String scriptsFolderName) {
		this.scriptsFolderName = scriptsFolderName;
	}

	public String getLogsFolder() {
		return logsFolder;
	}

	public void setLogsFolder(String logsFolder) {
		this.logsFolder = logsFolder;
	}

	public String getErrorsFolder() {
		return errorsFolder;
	}

	public void setErrorsFolder(String errorsFolder) {
		this.errorsFolder = errorsFolder;
	}

	public String getMapperComment() {
		return mapperComment;
	}

	public void setMapperComment(String mapperComment) {
		this.mapperComment = mapperComment;
	}

	public String getMapper() {
		return mapper;
	}

	public void setMapper(String mapper) {
		this.mapper = mapper;
	}

	public String getBatchSizeComment() {
		return batchSizeComment;
	}

	public void setBatchSizeComment(String batchSizeComment) {
		this.batchSizeComment = batchSizeComment;
	}

	public int getBatchSize() {
		return batchSize;
	}

	public void setBatchSize(int batchSize) {
		this.batchSize = batchSize;
	}

	public M getMapperVars() {
		return mapperVars;
	}

	public void setMapperVars(M mapperVars) {
		this.mapperVars = mapperVars;
	}

	public String getThreadsComment() {
		return threadsComment;
	}

	public void setThreadsComment(String threadsComment) {
		this.threadsComment = threadsComment;
	}

	public int getThreads() {
		return threads;
	}

	public void setThreads(int threads) {
		this.threads = threads;
	}

	public String getWalltimeComment() {
		return walltimeComment;
	}

	public void setWalltimeComment(String walltimeComment) {
		this.walltimeComment = walltimeComment;
	}

	public String getWalltime() {
		return walltime;
	}

	public void setWalltime(String walltime) {
		this.walltime = walltime;
	}

	public String getMaxMemoryComment() {
		return maxMemoryComment;
	}

	public void setMaxMemoryComment(String maxMemoryComment) {
		this.maxMemoryComment = maxMemoryComment;
	}

	public String getMaxMemory() {
		return maxMemory;
	}

	public void setMaxMemory(String maxMemory) {
		this.maxMemory = maxMemory;
	}

	public String getSlurmUserNameComment() {
		return slurmUserNameComment;
	}

	public void setSlurmUserNameComment(String slurmUserNameComment) {
		this.slurmUserNameComment = slurmUserNameComment;
	}

	public String getSlurmUserName() {
		return slurmUserName;
	}

	public void setSlurmUserName(String slurmUserName) {
		this.slurmUserName = slurmUserName;
	}

	public String getFinishedemailaddressComment() {
		return finishedemailaddressComment;
	}

	public void setFinishedemailaddressComment(String finishedemailaddressComment) {
		this.finishedemailaddressComment = finishedemailaddressComment;
	}

	public String getFinishedemailaddress() {
		return finishedemailaddress;
	}

	public void setFinishedemailaddress(String finishedemailaddress) {
		this.finishedemailaddress = finishedemailaddress;
	}

	public String getPairStringsComment() {
		return pairStringsComment;
	}

	public void setPairStringsComment(String pairStringsComment) {
		this.pairStringsComment = pairStringsComment;
	}

	public String[] getPairStrings() {
		return pairStrings;
	}

	public void setPairStrings(String[] pairStrings) {
		this.pairStrings = pairStrings;
	}
	@Override
	public HashMap<String, Integer> getStringAllowence() {//function to be overwritten by child classes
		HashMap<String, Integer> stringAllowence = super.getStringAllowence();
		stringAllowence.put("outModeComment", 0);
		stringAllowence.put("outMode", 0);
		stringAllowence.put("saveGenomeComment", 0);
		stringAllowence.put("saveGenome", 0);
			
		return stringAllowence;
	}
}
