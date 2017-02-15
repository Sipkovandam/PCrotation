package STAR;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.sun.msv.datatype.xsd.datetime.Util;

import Kallisto.Slurm;
import Tools.FileUtils;
import Tools.JSONutil;
import Tools.Runnable;
import Tools.Script;
import Tools.FileSearcher;

public class _STAR_Pipeline extends Script<_STAR_Pipeline> {
	
	// This class does the following:
	// 1. Finds all fastQ files in a folder and subfolders
	// 2. Determines which fastQ files are pairs
	// 3. Maps the files using Kallisto
	// 4. Creates a file describing the number of reads per sample mapped
	// 5. Sums reads of transcritps to genes
	// 6. Creates 1 large matrix for genes and one for transcripts
	// 7. Creates 1 large matrix for samples >70% mapping

	/**
	 * 
	 */
	private static final long serialVersionUID = 1053814892719523118L;
	// groups/umcg-wijmenga/tmp04/umcg-svandam/FeatureCounts/subread-1.5.1-source/bin/featureCounts
	private String comments = "//all variables in this file containing (Comment) are comments on the variable below and do not need to be initiated";
	private String passComment = "//1 for first pass only. 2 for first and second pass. -2 for secondpass only.";
	private int pass = 2;// 1 for first pass only. 2 for first and second pass.
							// -2 for secondpass only.
	private String outputFolderComment = "/root/directory/; MANDATORY //Folder where all output will be written (empty folder recommended)";
	private String outputFolder = null;// "E:/Groningen/Test/STAR/";

	FileSearcher fileSearcher = new FileSearcher(outputFolder);
	private String fastQfilesFNComment = "/root/directory/fileNamesFile.txt; OPTIONAL //Paths to file containing, enter separated absolute path containing, fastQ filenames to include in the analysis. OVERWRITES/REPLACES (inputFolder), (fastQSearchStrings) and (forbiddenSearchStrings)";
	private String fastQfilesFN = null;// list of paths to fastQ files to include in the analysis (overwrites inputFolder, fastQSearchStrings and forbiddenSearchStrings)

	// Mapper variables in case of star need to run: STAR --runMode genomeGenerate --genomeDir /groups/umcg-gdio/tmp04/umcg-svandam/Data/RNAseq/Annotation/GCRh38/STARindex/ --genomeFastaFiles /groups/umcg-gdio/tmp04/umcg-svandam/Data/RNAseq/Annotation/GCRh38/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa --runThreadN 12
	Slurm<STAR_Variables> slurmVars = new Slurm<STAR_Variables>();
	// SpliceSites variables
	SpliceMerger spliceMerger = new SpliceMerger();

	@Override
	public void run() {
		try {
			List<Runnable> steps = initiate();

			// run the selected steps
			for (Runnable step : steps)
				step.run();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private List<Runnable> initiate() throws CloneNotSupportedException {
		mkdir(getOutputFolder());

		// set slurm variable used in each run
		slurmVars.setOutputRoot(this.getOutputFolder());

		// set fileSearcher arguments
		fileSearcher.setWriteName(this.getOutputFolder() + "fastQfiles.txt");
		fileSearcher.setJsonFN(this.getOutputFolder()+fileSearcher.getNewJsonFN());
		fileSearcher.writeConfig();
		final String fastQFiles = getFastq();
		this.slurmVars.setFastQFiles(fastQFiles);
		
		
		// First pass genome build arguments (adding in new splice sites)
		String firstPassGenomeDir = getOutputFolder() + "firstPassGenome/";
		Slurm<STAR_Variables> firstBuildVars = (Slurm<STAR_Variables>) slurmVars.clone(this.getOutputFolder() + this.slurmVars.getNewJsonFN("firstPassGenomeBuild"));
		firstBuildVars.setOutputRoot(firstPassGenomeDir);
		firstBuildVars.mapperVars.setMaxMemoryGenomeBuild(this.slurmVars.mapperVars.getMaxMemoryGenomeBuild());
		firstBuildVars.mapperVars.setSTAR_Extra_Arguments(
				slurmVars.mapperVars.getSTAR_Extra_Arguments() + " --sjdbInsertSave All --limitSjdbInsertNsj 10000000");
		firstBuildVars.mapperVars.setSjdbFileChrStartEnd(slurmVars.mapperVars.getSjdbFileChrStartEnd());
		firstBuildVars.mapperVars.setOutMode("None");
		firstBuildVars.mapperVars.setSaveGenome(true);
		firstBuildVars.mapperVars.setCountExpression(true);// might come in handy to
		firstBuildVars.writeConfig();

		// first pass STAR arguments
		Slurm<STAR_Variables> sTARv_FirstPass = (Slurm<STAR_Variables>) slurmVars
				.clone(this.getOutputFolder() + this.slurmVars.getNewJsonFN("FirstPass"));
		sTARv_FirstPass.mapperVars.setGenomeDir(firstPassGenomeDir+"Results/_STARgenome/");
		sTARv_FirstPass.setOutputRoot(this.slurmVars.getOutputRoot() + "1stPass/");
		sTARv_FirstPass.set_STAR_Folder(sTARv_FirstPass.getOutputRoot() + "Results/");
		sTARv_FirstPass.mapperVars.setOutMode("None");
		sTARv_FirstPass.mapperVars.setSjdbFileChrStartEnd(null);
		sTARv_FirstPass.mapperVars.setgTFfile(null);
		sTARv_FirstPass.writeConfig();
		// first pass merge splice junction files arguments
		SpliceMerger mergeV_FirstPass = (SpliceMerger) this.spliceMerger
				.clone(this.getOutputFolder() + this.spliceMerger.getNewJsonFN("FirstPass"));
		mergeV_FirstPass.setWriteFN_Splice(
				FileUtils.makeFolderNameEndWithSlash(this.getOutputFolder()) + "SJ_Merged_1stPass.out.tab");
		mergeV_FirstPass.setInputFolder_Splice(sTARv_FirstPass.get_STAR_Folder());
		mergeV_FirstPass.setWriteFolder_SplicePerGene(sTARv_FirstPass.getOutputRoot());
		mergeV_FirstPass.writeConfig();

		// Second pass rebuild genome arguments (adding in new splice sites)
		String secondPassGenomeDir = getOutputFolder() + "secondPassGenome/";
		Slurm<STAR_Variables> secondBuildVars = (Slurm<STAR_Variables>) slurmVars.clone(this.getOutputFolder() + this.slurmVars.getNewJsonFN("secondPassGenomeBuild"));
		secondBuildVars.setOutputRoot(secondPassGenomeDir);
		secondBuildVars.mapperVars.setMaxMemoryGenomeBuild(this.slurmVars.mapperVars.getMaxMemoryGenomeBuild());
		secondBuildVars.mapperVars.setSTAR_Extra_Arguments(
				slurmVars.mapperVars.getSTAR_Extra_Arguments() + " --sjdbInsertSave All --limitSjdbInsertNsj 10000000");
		secondBuildVars.mapperVars.setSjdbFileChrStartEnd(mergeV_FirstPass.getWriteFN_Splice());
		secondBuildVars.mapperVars.setOutMode("None");
		secondBuildVars.mapperVars.setSaveGenome(true);
		secondBuildVars.mapperVars.setCountExpression(true);// might come in handy to
		secondBuildVars.writeConfig();												// check things

		// second pass STAR arguments
		Slurm<STAR_Variables> sTARv_SecondPass = (Slurm<STAR_Variables>) slurmVars.clone(this.getOutputFolder() + this.slurmVars.getNewJsonFN("SecondPass"));
		sTARv_SecondPass.mapperVars.setGenomeDir(secondPassGenomeDir+"Results/_STARgenome/");
		sTARv_SecondPass.setOutputRoot(this.slurmVars.getOutputRoot() + "2ndPass/");
		sTARv_SecondPass.set_STAR_Folder(sTARv_SecondPass.getOutputRoot() + "Results/");
		sTARv_SecondPass.mapperVars.setOutMode("BAM Unsorted");
		sTARv_SecondPass.mapperVars.setCountExpression(true);
		sTARv_SecondPass.mapperVars.setSjdbFileChrStartEnd(null);
		sTARv_SecondPass.writeConfig();
		// second pass merge splice junction files arguments
		SpliceMerger mergeV_SecondPass = (SpliceMerger) this.spliceMerger
				.clone(this.getOutputFolder() + this.spliceMerger.getNewJsonFN("SecondPass"));
		mergeV_SecondPass.setWriteFN_Splice(sTARv_SecondPass.getOutputRoot() + "SJ_Merged_2ndPass.out.tab");
		mergeV_SecondPass.setInputFolder_Splice(sTARv_SecondPass.get_STAR_Folder());
		mergeV_SecondPass.setWriteFolder_SplicePerGene(sTARv_SecondPass.getOutputRoot());
		mergeV_SecondPass.writeConfig();

		// create exon expression files and calculate expression ratios per gene
		ExonExpressionMerger exonExpressionMerger = new ExonExpressionMerger();
		exonExpressionMerger.setJsonFN(this.getOutputFolder() + exonExpressionMerger.getNewJsonFN("SecondPass"));
		exonExpressionMerger.setSearchFolder(sTARv_SecondPass.getOutputRoot());
		exonExpressionMerger.setWriteFolder(sTARv_SecondPass.getOutputRoot());
		exonExpressionMerger.setEnsgToGeneSymbolFn(spliceMerger.spliceRatioCalculator.getEnsgToGeneSymbolFN());
		exonExpressionMerger.writeConfig();
		
		// add steps (in right order of course ;))
		List<Runnable> steps = new ArrayList<>();
		if (this.getFastQfilesFN() == null)// if list of fastq files is not
											// supplied yet, find them in the
											// supplied folder(s)
			steps.add(this.fileSearcher);
		if (this.getPass() == 1 || this.getPass() == 2) {
			steps.add(firstBuildVars);
			steps.add(sTARv_FirstPass);
			steps.add(mergeV_FirstPass);
		}
		if (this.getPass() == 2 || this.getPass() == -2) {
			steps.add(secondBuildVars);
			steps.add(sTARv_SecondPass);
			steps.add(mergeV_SecondPass);
			steps.add(exonExpressionMerger);
		}
		return steps;

	}

	private void mkdir(String folderName) {
		if (folderName == null) {
			p("folderName is null: please set variable folderName\t");
			System.exit(2);
		}
		File folder = new File(folderName);
		File parent = new File(folder.getParent());
		if (!parent.exists()) {
			p("Exiting, parent folder" + parent.getAbsolutePath() + " does not exist:\t");
			System.exit(2);
		}

		if (!new File(getOutputFolder()).exists())
			folder.mkdir();
	}

	private String getFastq() {
		if (this.getFastQfilesFN() != null)
			return this.getFastQfilesFN();;
		return fileSearcher.getWriteName();
	}

	public String includeWritePath(String name) {
		if (!name.contains("\\") && !name.contains("/"))
			return FileUtils.makeFolderNameEndWithSlash(this.outputFolder) + name;
		return name;
	}

	public static _STAR_Pipeline readVars(String jsonFilename, String executeClassName)// function
																						// should
																						// be
																						// in
																						// separate
																						// interface...
	{
		String jsonString = "";
		try {
			File file = new File(jsonFilename);
			if (!file.exists()) {
				System.out.println("Json file does not exist:" + jsonFilename);
				System.exit(0);
			}
			jsonString = new String(Files.readAllBytes(Paths.get(jsonFilename)));
		} catch (Exception e) {
			e.printStackTrace();
		}

		Gson gson = new GsonBuilder().setPrettyPrinting().create();
		_STAR_Pipeline var = gson.fromJson(jsonString, _STAR_Pipeline.class);
		var.jsonFN = jsonFilename;
		System.out.println(gson.toJson(var));
		return var;
	}

	public String getOutputFolder() {
		String outputFolder = FileUtils.makeFolderNameEndWithSlash(this.outputFolder);
		return outputFolder;
	}

	public _STAR_Pipeline clone() throws CloneNotSupportedException {
		_STAR_Pipeline copy = (_STAR_Pipeline) super.clone();
		return copy;
	}

	public int getPass() {
		return pass;
	}

	private void setPass(int pass) {
		this.pass = pass;
	}

	public String getFastQfilesFN() {
		return fastQfilesFN;
	}

	private void setFastQfilesFN(String fastQfilesFN) {
		this.fastQfilesFN = fastQfilesFN;
	}

	public void createConfig(String jsonFN) {
		this.slurmVars.mapperVars = new STAR_Variables();
		createConfig(jsonFN, this, false, true);
	}
	@Override
	public HashMap<String, Integer> getStringAllowence() {//function to be overwritten by child classes
		HashMap<String, Integer> stringAllowence = super.getStringAllowence();
		stringAllowence.put("writeNameComment", 0);
		stringAllowence.put("writeName", 0);
		stringAllowence.put("fastQFiles", 0);
		stringAllowence.put("outputRoot", 0);
		stringAllowence.putAll(this.slurmVars.getStringAllowence());
		stringAllowence.put("inputFolder_SpliceComment", 0);
		stringAllowence.put("inputFolder_Splice", 0);
		stringAllowence.put("writeFN_SpliceComment", 0);
		stringAllowence.put("writeFN_Splice", 0);
		stringAllowence.put("writeFolder_SplicePerGeneComment", 0);
		stringAllowence.put("writeFolder_SplicePerGene", 0);
		
		return stringAllowence;
	}
	
	
}