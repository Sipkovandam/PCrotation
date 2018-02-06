package Slurm;

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

import Kallisto.Kallisto_ClusterHandler;
import STAR.FastqDeduplicator;
import STAR.FastqDeduplicator_ClusterHandler;
import STAR.STAR_ClusterHandler;
import Tests.ReadIncludedFile;
import Tools.ExecCommand;
import Tools.FileUtils;
import Tools.JSONutil;
import Tools.Script;
import Tools.Util;

public class ClusterHandler <M> extends Script<ClusterHandler<M>> implements Cloneable{

	//@Rule
	//public TemporaryFolder tempFolder = new TemporaryFolder();

	/**
	 * 
	 */	
	private static final long serialVersionUID = 1L;
	private String fastQFiles = null;
	private String outputRoot = null;
	private transient String scriptsFolderName = null;
	private transient String logsFolder = null;
	private transient String errorsFolder = null;
	private transient String sTAR_Folder = null;
	private String jobNameComment = "STAR/2.5.1b-foss-2015b or Kallisto/0.42.2.1-goolf-1.7.20;MANDATORY // The executeable of the STAR version to be used";
	private String jobName = "STAR/2.5.1b-foss-2015b";// "STAR/2.5.1b-foss-2015b";
	private String batchSizeComment = "MANDATORY //number of files to run in 1 batch. Each batch only loads the genome once to save runtime. However if the number is to large the walltime may be hit and some files will not/partly be included in the analysis; no warnings are generated when this happens but if the walltime has been hit can be checked in the /error folder";
	private int batchSize = 1;
	//mapper variables
	public M jobType = null;
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
	
	private String clusterHandlerComment = "slurm for slurm, shark for shark";
	private String clusterHandler = "slurm";

	// pairStrings to be in pairs. String in uneven index is replaced with that
	// in the even index to get the paired end partner.
	// All file names containing a string in one of the uneven indexes is
	// skipped)
	
	public ClusterHandler(String jobName)
	{
		this.setMapper(jobName);
		if(isSTAR())
			jobType = (M) new STAR_ClusterHandler();
		else if(isKallisto())
			jobType = (M) new Kallisto_ClusterHandler();
		else if(isFastqDeduplicator())
			jobType = (M) new FastqDeduplicator_ClusterHandler();
		else
		{
			log("Please select an appropriate mapper containing \"FastqDeduplicator\",\"STAR\" of \"kallisto\" in the jobName name");
		}
	}

	public void run() {
		try
		{
			//checkArgs(args);
			init();
			((SlurmJob) jobType).createSlurmFiles(fastQFiles, this);
	
			runSlurmScripts();
		}catch(Exception e){e.printStackTrace();}
	}

	private void init() {
		if (outputRoot == null && fastQFiles!=null)
			outputRoot = new File(fastQFiles).getParent() + "/";
		log("WriteFolder=" +outputRoot);
		
		scriptsFolderName = outputRoot + "Scripts/";
		logsFolder = outputRoot + "Logs/";
		errorsFolder = outputRoot + "Errors/";
		if(sTAR_Folder==null && outputRoot!=null)
			sTAR_Folder = outputRoot + "Results/";
		
		FileUtils.makeDir(scriptsFolderName);
		FileUtils.makeDir(logsFolder);
		FileUtils.makeDir(errorsFolder);
		FileUtils.makeDir(sTAR_Folder);
		
		this.writeConfig();
	}
	
	private boolean isFastqDeduplicator()
	{
		return jobName.toLowerCase().startsWith("bbmap");
	}
	
	private boolean isSTAR() {
		return jobName.toLowerCase().contains("star");
	}
	
	private boolean isKallisto() {
		return jobName.toLowerCase().contains("kallisto");
	}

	public void writeConfig(String jsonFN)
	{
		String kallistoFN = FileUtils.addBeforeExtention(jsonFN, "_kallisto");
		this.jobType= (M) new Kallisto_ClusterHandler();
		super.writeConfig(kallistoFN, this);
		String STAR_FN = FileUtils.addBeforeExtention(jsonFN, "_STAR");
		this.jobType= (M) new STAR_ClusterHandler();
		super.writeConfig(STAR_FN, this);
	}
	
	private void runSlurmScripts() throws FileNotFoundException, IOException {
		// carefull with this. Kallisto200sh script does not update! You need to
		// delete and reimport project of that!

		String shellFN  = getShellFileToRun();
		String scriptTrackingFile = createScriptsTrackingFile();
		
		String mappingPercentagesFN = new File(sTAR_Folder).getParent() + "/mappingPerSample.txt";
		String command = "bash " + shellFN + " " + scriptsFolderName + "*.sh " + sTAR_Folder + " "
				+ slurmUserName + " " + mappingPercentagesFN + " " + finishedemailaddress + " " + jobName + " " + scriptTrackingFile;
		
		log("Output folder of slurm scripts:\t" + this.getScriptsFolderName());
		log("Running SLURM scripts using:\t" + command);	
		
//		System.exit(2);
		
		ExecCommand exec = new ExecCommand(command);
		exec.reportErrors();
	}
	private String getShellFileToRun()
	{
		String shellFN = null; 
		if(isSlurmCluster())
			shellFN=getSlurmShellFile();
		else if(isSharkCluster())
			shellFN=getSharkShellFile();
			
		FileUtils.extractFile(ClusterHandler.class.getProtectionDomain().getCodeSource().getLocation().getPath(), shellFN);
		
		if(!new File(shellFN).exists())
		{
			log("Slurmscript does not exist, please copy to:\t" + shellFN);
			System.exit(2);
		}
		return shellFN;
	}

	private String getSharkShellFile()
	{
		String shellFN = "resources/Kallisto_Shark.sh";//FileUtils.prepareBinaryFromJar("Kallisto200.sh");
		if(isSTAR())
			shellFN="resources/STAR_OnePerNode_MultiUser_Shark.sh";
		return shellFN;
	}

	private String getSlurmShellFile()
	{
		String shellFN = "resources/Kallisto_Slurm.sh";//FileUtils.prepareBinaryFromJar("Kallisto200.sh");
		if(isSTAR())
			shellFN="resources/STAR_OnePerNode_MultiUser_Slurm.sh";
		return shellFN;
	}

	private String createScriptsTrackingFile()
	{
		String startedScriptsFn = "scriptsStarted.txt";
		//remove any old tracking file.
		new File(startedScriptsFn).delete();
		return startedScriptsFn;
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


	private STAR_ClusterHandler getSTAR_V() {
		return (STAR_ClusterHandler) jobType;
	}
	public String getOutputRoot() {
		return FileUtils.makeFolderNameEndWithSlash(outputRoot);
	}

	public void setOutputRoot(String outputRoot) {
		this.outputRoot = outputRoot;
	}

	private Kallisto_ClusterHandler getKallistoV() {
		return (Kallisto_ClusterHandler) jobType;
	}


	public ClusterVariables writeSlurmCommands(BufferedWriter writer, int fileNumber) throws Exception {
		ClusterVariables clusterCommands = null;
		String jN = getJobName(jobName);
		if(isSlurmCluster())
			clusterCommands = new SlurmVariables();
		else if(isSharkCluster())
			clusterCommands = new SharkVariables();
		
		writeClusterCommand(writer,
							true,
							clusterCommands.getFileHeader() + "\n");
		writeClusterCommand(writer,
							true,
							clusterCommands.getJobName() + "j" + jN + fileNumber + ".sh" + "\n");
		writeClusterCommand(writer,
							true,
							clusterCommands.getLogsFolder() + logsFolder + fileNumber + ".out\n");
		writeClusterCommand(writer,
							true,
							clusterCommands.getErrorsFolder() + errorsFolder + fileNumber + ".err\n");
		writeClusterCommand(writer,
							true,
							clusterCommands.getWalltime() + walltime + "\n");
		writeClusterCommand(writer,
							true,
							clusterCommands.getThreads() + threads + "\n");
		writeClusterCommand(writer,
							true,
							clusterCommands.getMaxMemory() + maxMemory + "\n");

		writeClusterCommand(writer,
							true,
							clusterCommands.getExtra());
		writeClusterCommand(writer,
							true,
							clusterCommands.getLoadModule() + jobName + clusterCommands.getLoadModule2() + "\n");

		return clusterCommands;
	}
	
	private void writeClusterCommand(	BufferedWriter writer,
	                                 	boolean conditionPass,
										String command) throws IOException
	{
		if(conditionPass)
			writer.write(command+ "\n");
		
	}

	public boolean isSlurmCluster()
	{
		return clusterHandler.toLowerCase().equals("slurm");
	}

	public boolean isSharkCluster()
	{
		return clusterHandler.toLowerCase().equals("shark");
	}

	private String getJobName(String jobName)
	{
		if(clusterHandler.toLowerCase().equals("slurm"))
			return jobName;
		else if(clusterHandler.toLowerCase().equals("shark"))
		{
			if(jobName.toLowerCase().contains("star"))
				return "starJob";
			if(jobName.toLowerCase().contains("kallisto"))
				return "kallistoJob";
		}
		return jobName;
	}

	void checkArgs(String[] args) {
		if (System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko") && args.length < 1)
			return;
		if (args.length < 1) {
			log("Script requires the following arguments:\n"
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
			log("arg=" + args[a]);
			String[] split =args[a].split("=");
			String arg = split[0];
			if(split.length <2)
			{
				log("Argument missing, skipping:" + args[a]);
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
				jobName = parseString(value);
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
				jobName = parseString(value);
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
				log("Incorrect argument supplied:\n" + args[a] + "\nexiting");
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
		return jobNameComment;
	}

	public void setMapperComment(String mapperComment) {
		this.jobNameComment = mapperComment;
	}

	public String getMapper() {
		return jobName;
	}

	public void setMapper(String mapper) {
		this.jobName = mapper;
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
		return jobType;
	}

	public void setMapperVars(M mapperVars) {
		this.jobType = mapperVars;
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

	@Override
	public HashMap<String, Integer> getStringAllowence() {//function to be overwritten by child classes
		HashMap<String, Integer> stringAllowence = super.getStringAllowence();
		stringAllowence.put("outModeComment", 0);
		stringAllowence.put("outMode", 0);
		stringAllowence.put("saveGenomeComment", 0);
		stringAllowence.put("saveGenome", 0);
			
		return stringAllowence;
	}


	public String getSTAR_Folder()
	{
		return sTAR_Folder;
	}


	public void setSTAR_Folder(String sTAR_Folder)
	{
		this.sTAR_Folder = sTAR_Folder;
	}

	public String getClusterHandler()
	{
		return clusterHandler;
	}

	public void setClusterHandler(String clusterHandler)
	{
		this.clusterHandler = clusterHandler;
	}
}
