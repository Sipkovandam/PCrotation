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

import MatrixScripts.RowStatisticsGetter;
import MatrixScripts.CenterPerRow;
import MatrixScripts.LaneMerger;
import MatrixScripts.LargeFileSorter;
import MatrixScripts.MyMatrix;
import MatrixScripts.Transpose;
import PCA.CorrelationLarge;
import PCA.Pca;
import PCA.PcaPipelineLite;
import Slurm.ClusterHandler;
import TextEditing.FileSplitter;
import Tools.FileUtils;
import Tools.JSONutil;
import Tools.Runnable;
import Tools.Script;
import Tools.FileSearcher;

public class STAR_Pipeline extends Script<STAR_Pipeline>
{

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

	FileSearcher fileSearcherFirstPass = new FileSearcher(outputFolder);
	private String fastQfilesFirstPassFNComment = "/root/directory/fileNamesFile.txt; OPTIONAL //Paths to file containing, enter separated absolute path containing, fastQ filenames to include in the first pass of the analysis. OVERWRITES/REPLACES (inputFolder), (fastQSearchStrings) and (forbiddenSearchStrings)";
	private String fastQfilesFirstPassFN = null;// list of paths to fastQ files to include in the analysis (overwrites inputFolder, fastQSearchStrings and forbiddenSearchStrings)

	FileSearcher fileSearcherSecondPass = new FileSearcher(outputFolder);
	private String fastQfilesSecondPassFNComment = "/root/directory/fileNamesFile.txt; OPTIONAL //Paths to file containing, enter separated absolute path containing, fastQ filenames to include in the second pass of the analysis. OVERWRITES/REPLACES (inputFolder), (fastQSearchStrings) and (forbiddenSearchStrings)";
	private String fastQfilesSecondPassFN = null;// list of paths to fastQ files to include in the analysis (overwrites inputFolder, fastQSearchStrings and forbiddenSearchStrings)

	private String sampleSheetFnComment = "/root/directory/samplesheet.txt; OPTIONAL;//";
	private String sampleSheetFn = null;
	// Mapper variables in case of star need to run: STAR --runMode genomeGenerate --genomeDir /groups/umcg-gdio/tmp04/umcg-svandam/Data/RNAseq/Annotation/GCRh38/STARindex/ --genomeFastaFiles /groups/umcg-gdio/tmp04/umcg-svandam/Data/RNAseq/Annotation/GCRh38/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa --runThreadN 12
	ClusterHandler<STAR_ClusterHandler> slurmVars = new ClusterHandler<STAR_ClusterHandler>("STAR/2.5.1b-foss-2015b");
	//GeneNameAdder variables
	private String annotationFNComment = "/root/directory/annotation.txt; MANDATORY // A file containing ensembl IDs in the 1st column, chromosome in 2nd, start position in 3rd, end position in 4th column";
	private String annotationFN = null;
	private String ensgToGeneSymbolFNComment = "/root/directory/ensemblToGeneSymbolFN.txt; MANDATORY // A file containing ensembl IDs in the first column and corresponding Gene Symbols in the second column";
	private String ensgToGeneSymbolFN = null;//"E:/Groningen/Data/Annotation/GRCh38/EnsgToGeneSymbol.txt";//"/groups/umcg-gdio/tmp04/umcg-svandam/Data/RNAseq/Annotation/GCRh38/EnsgToGeneSymbol.txt";//


	/** Runs the pipeline entailing the following steps
	 *  1. Find all the fastQ files that should be included in the first step
	 *  2. Build the first pass STAR genome annotation including all defined splice junctions
	 *  3. Run the first pass mapping on all the samples to be included (from step 1)
	 *  4. Merge the splice files and select the junctions to be used in the second pass STAR run
	 *  5. Find all the fastQ files that should be included in the second step
	 *  6. Build the second pass STAR genome annotation including all junctions selected from the first run
	 *  7. Run the second pass mapping on all the samples to be included (from step 5)
	 *  8. Merge the splice files and select the junctions detected in the second run. Also create a file per gene containing all splice junction counts for all samples for that gene.
	 *  9. Merge the exon files.
	 */

	@Override
	public void run()
	{
		try
		{
			List<Runnable> steps = initiate();

			// run the selected steps
			MyMatrix passMatrix = null;
			for (Runnable step : steps)
				passMatrix=step.run(passMatrix);
		} catch (Exception e)
		{
			e.printStackTrace();
		}
	}

	private List<Runnable> initiate() throws CloneNotSupportedException
	{
		mkdir(getOutputFolder());

		// set slurm variable used in each run
		slurmVars.setOutputRoot(this.getOutputFolder());

		// set fileSearcher arguments
		fileSearcherFirstPass.setWriteName(this.getOutputFolder() + "fastQfilesFirstPass.txt");
		fileSearcherFirstPass.setJsonFN(this.getOutputFolder() + fileSearcherFirstPass.getNewJsonFN());
		fileSearcherFirstPass.writeConfig();
		
		final String fastQFilesFirstPass = getFastqFirstPass();
		this.slurmVars.setFastQFiles(fastQFilesFirstPass);

		// First pass genome build arguments (adding in new splice sites)
		String firstPassGenomeDir = getOutputFolder() + "firstPassGenome/";
		ClusterHandler<STAR_ClusterHandler> firstBuildVars = (ClusterHandler<STAR_ClusterHandler>) slurmVars.clone(this.getOutputFolder() + this.slurmVars.getNewJsonFN("firstPassGenomeBuild"));
		firstBuildVars.setOutputRoot(firstPassGenomeDir);
		firstBuildVars.jobType.setMaxMemoryGenomeBuild(this.slurmVars.jobType.getMaxMemoryGenomeBuild());
		firstBuildVars.jobType.setSTAR_Extra_Arguments(slurmVars.jobType.getSTAR_Extra_Arguments() + " --sjdbInsertSave All --limitSjdbInsertNsj 10000000");
		firstBuildVars.jobType.setSjdbFileChrStartEnd(slurmVars.jobType.getSjdbFileChrStartEnd());
		firstBuildVars.jobType.setOutMode("None", null);
		firstBuildVars.jobType.setSaveGenome(true);
		firstBuildVars.jobType.setCountExpression(true);// might come in handy to
		firstBuildVars.writeConfig();

		// first pass STAR arguments
		ClusterHandler<STAR_ClusterHandler> sTARv_FirstPass = (ClusterHandler<STAR_ClusterHandler>) slurmVars.clone(this.getOutputFolder() + this.slurmVars.getNewJsonFN("FirstPass"));
		sTARv_FirstPass.jobType.setGenomeDir(firstPassGenomeDir + "Results/_STARgenome/");
		sTARv_FirstPass.setOutputRoot(this.slurmVars.getOutputRoot() + "1stPass/");
		sTARv_FirstPass.set_STAR_Folder(sTARv_FirstPass.getOutputRoot() + "Results/");
		sTARv_FirstPass.jobType.setOutMode("None", null);
		sTARv_FirstPass.jobType.setSjdbFileChrStartEnd(null);
		sTARv_FirstPass.jobType.setgTFfile(null);
		sTARv_FirstPass.getMapperVars().setPicardJarFn(null);
		sTARv_FirstPass.getMapperVars().setSamtoolsVersion(null);
		sTARv_FirstPass.getMapperVars().setOutMode("None", null);
		sTARv_FirstPass.writeConfig();
		
		//merge files ran on multiple lanes
		LaneMerger_StarOutput laneMerger_StarOuptut_FirstPass = new LaneMerger_StarOutput(); 
		laneMerger_StarOuptut_FirstPass.setStarResultsFolderName(sTARv_FirstPass.get_STAR_Folder());
		
		// first pass merge splice junction files arguments
		SpliceMerger mergeV_FirstPass = new SpliceMerger();
		mergeV_FirstPass.setspliceFilesFolder(sTARv_FirstPass.get_STAR_Folder());
		mergeV_FirstPass.setSpliceFileExtencion("SJ.out.tab.gz");
		mergeV_FirstPass.setAnnotatedOnly(false);
		mergeV_FirstPass.setWriteFn(sTARv_FirstPass.getOutputRoot() + "Expression_SJ_Merged_1stPass.txt.gz");
		mergeV_FirstPass.setReadSumCutoff(8);
		mergeV_FirstPass.setExonsInsteadOfJunctions(false);
		mergeV_FirstPass.setWriteFn_Splice_2ndPassInput(sTARv_FirstPass.getOutputRoot() + "Expression_SJ_Merged_1stPass_readCutoff"+mergeV_FirstPass.getReadSumCutoff()+".txt");
		mergeV_FirstPass.setWriteSpliceSummaryFn(sTARv_FirstPass.getOutputRoot() + "Expression_SJ_Summary_1stPass.txt");

		// set fileSearcher arguments for the second pass
		fileSearcherSecondPass.setWriteName(this.getOutputFolder() + "fastQfilesSecondPass.txt");
		fileSearcherSecondPass.setJsonFN(this.getOutputFolder() + fileSearcherSecondPass.getNewJsonFN());
		fileSearcherSecondPass.writeConfig();
		final String fastQFilesSecondPass = getFastqSecondPass();
		this.slurmVars.setFastQFiles(fastQFilesSecondPass);

		// Second pass rebuild genome arguments (adding in new splice sites)
		String secondPassGenomeDir = getOutputFolder() + "secondPassGenome/";
		ClusterHandler<STAR_ClusterHandler> secondBuildVars = (ClusterHandler<STAR_ClusterHandler>) slurmVars.clone(this.getOutputFolder() + this.slurmVars.getNewJsonFN("secondPassGenomeBuild"));
		secondBuildVars.setOutputRoot(secondPassGenomeDir);
		secondBuildVars.jobType.setMaxMemoryGenomeBuild(this.slurmVars.jobType.getMaxMemoryGenomeBuild());
		secondBuildVars.jobType.setSTAR_Extra_Arguments(slurmVars.jobType.getSTAR_Extra_Arguments() + " --sjdbInsertSave All --limitSjdbInsertNsj 10000000");
		secondBuildVars.jobType.setSjdbFileChrStartEnd(mergeV_FirstPass.getWriteFN_Splice_2ndPassInput());//all variants used in first pass + all variants that are above the read cutoff minimum
		secondBuildVars.jobType.setOutMode("None", null);
		secondBuildVars.jobType.setSaveGenome(true);
		secondBuildVars.jobType.setCountExpression(true);// might come in handy to
		secondBuildVars.writeConfig(); // check things

		// second pass STAR arguments
		ClusterHandler<STAR_ClusterHandler> sTARv_SecondPass = (ClusterHandler<STAR_ClusterHandler>) slurmVars.clone(this.getOutputFolder() + this.slurmVars.getNewJsonFN("SecondPass"));
		sTARv_SecondPass.jobType.setGenomeDir(secondPassGenomeDir + "Results/_STARgenome/");
		sTARv_SecondPass.setOutputRoot(this.slurmVars.getOutputRoot() + "2ndPass/");
		sTARv_SecondPass.set_STAR_Folder(sTARv_SecondPass.getOutputRoot() + "Results/");
		String deduplicateAddString =  "deduplicated";
		sTARv_SecondPass.jobType.setOutMode("BAM SortedByCoordinate", deduplicateAddString);
		sTARv_SecondPass.jobType.setsTAR_Extra_Arguments("--limitBAMsortRAM 30000000000");
		sTARv_SecondPass.jobType.setCountExpression(true);
		sTARv_SecondPass.jobType.setSjdbFileChrStartEnd(null);
		sTARv_SecondPass.writeConfig();
		
		//merge files ran on multiple lanes
		LaneMerger_StarOutput laneMerger_StarOuptut_SecondPass = new LaneMerger_StarOutput(); 
		laneMerger_StarOuptut_SecondPass.setStarResultsFolderName(sTARv_SecondPass.get_STAR_Folder());
		
		// second pass merge splice junction files arguments
		SpliceMerger mergeV_SecondPass = new SpliceMerger();
		mergeV_SecondPass.setspliceFilesFolder(sTARv_SecondPass.get_STAR_Folder());
		mergeV_SecondPass.setSpliceFileExtencion("SJ.out.tab.gz");
		mergeV_SecondPass.setColumnHeaderStringToUse(2);
		mergeV_SecondPass.setRequiredStrings(new String[]{deduplicateAddString});
		mergeV_SecondPass.setWriteFn(sTARv_SecondPass.getOutputRoot() + "Expression_SJ_Merged_2ndPass.txt.gz");
		
//		//laneMerger arguments; Merges columns in file that belong to the same sample but are ran on multiple lanes
//		LaneMerger laneMerger = new LaneMerger();
//		laneMerger.setLaneIndexInFileName(1);
//		laneMerger.setCountsFn(mergeV_SecondPass.getWriteFn());
//		laneMerger.setSampleSheetFn(sampleSheetFn);
//		laneMerger.setSummedFn(FileUtils.removeExtention(laneMerger.getCountsFn()) + "_LanesMerged.txt.gz");
		
		//GeneNameAdder; Adds gene names to the splice sites
		GeneNameAdder geneNameAdder = new GeneNameAdder();
		geneNameAdder.setGtfFn(slurmVars.jobType.getGTFfile());
		geneNameAdder.setSpliceFn(mergeV_SecondPass.getWriteFn());
		geneNameAdder.setWriteFn(FileUtils.addBeforeExtention(geneNameAdder.getSpliceFn(),"_geneNamesAdded"));
		geneNameAdder.setWriteFnNoGeneSymbol(FileUtils.removeExtention(geneNameAdder.getWriteFn())+ "_1rowNameColumn.txt");
					
		//FileSplitter;
		FileSplitter fileSplitter = new FileSplitter();
		fileSplitter.setFn(geneNameAdder.getWriteFn());
		fileSplitter.setGzipFiles(true);
		fileSplitter.setSpliceNamesCol(2);
		fileSplitter.setSplitRegex(",");
		fileSplitter.setSplitFirstColumnInstead(false);
		fileSplitter.setWriteFolder(new File(geneNameAdder.getWriteFn()).getParent()+"/PerGene/");
		
		//RatioCalculator; calculate relative usage of each splice site for each gene
		RatioCalculator ratioCalculator = new RatioCalculator();
		ratioCalculator.setExpressionFolder(fileSplitter.getWriteFolder());
		ratioCalculator.setExpressionSortedWriteFn(FileUtils.addBeforeExtention(geneNameAdder.getWriteFnNoGeneSymbol(), "_sorted"));
		ratioCalculator.setWriteFn(new File(geneNameAdder.getWriteFn()).getParent()+"/ratios.txt.gz");
		
		//Transpose; transpose the matrix so the slice variants are on the columns
		Transpose transpose = new Transpose();
		transpose.setFileName(ratioCalculator.getWriteFn());
		transpose.setWriteFn(FileUtils.removeExtention(transpose.getFileName())+ "_transposed.txt.gz");
		
		String pcaDir = FileUtils.removeExtention(transpose.getFileName())+"Pca/";
		new File(pcaDir).mkdirs();
		//CenterRows
		CenterPerRow centerPerRow = new CenterPerRow();
		centerPerRow.setFileName(transpose.getWriteFn());
		centerPerRow.setAveragesWriteFn(FileUtils.removeExtention(transpose.getWriteFn())+"_sampleAverages.txt");
		centerPerRow.setWriteFn(pcaDir+"MATRIX_Centered.txt.gz");
		
		//CorrelationLarge; calculate a covariance matrix over the samples dimension in the ratios file;
		CorrelationLarge correlationLarge = new CorrelationLarge();
		correlationLarge.setExpressionFN(centerPerRow.getWriteFn());
		correlationLarge.setCorrelation(false);
		correlationLarge.setThreads(slurmVars.getThreads());
		correlationLarge.setWriteFn(pcaDir+"ratios_rowCentered_covariance.txt.gz");
		
		//run PCA over the covariance matrix
		Pca pca = new Pca();
		pca.setInputMatrix(correlationLarge.getWriteFn());
		pca.setEigenVectorWriteFn(pcaDir+ "SAMPLE.eigenvectors.txt");
		pca.setEigenValueWriteFn(pcaDir+"SAMPLE.eigenvalues.txt");
		
		//calculate the averages per gene (not used for calculations, just to indicate which genes should be included)
		RowStatisticsGetter averagesPerRow = new RowStatisticsGetter();
		averagesPerRow.setFileName(ratioCalculator.getWriteFn());
		averagesPerRow.setWriteFn(pcaDir+"SAMPLE_Norm_GeneAverages.txt");
		averagesPerRow.setAbsolute(false);
		
		//run PC correction on the ratios file
		String correctedDir = pcaDir+"ratiosCorrected/";
		PcaPipelineLite PCApipelineLite = new PcaPipelineLite();
		PCApipelineLite.setSampleFile(ratioCalculator.getWriteFn());
		PCApipelineLite.setRunMode(2);
		PCApipelineLite.setPcaOverGenes(false);//not used, it uses centered matrix from PCA
		PCApipelineLite.setCenterGenes(false);//not used, it uses centered matrix from PCA
		PCApipelineLite.setCenterSamples(true);//not used, it uses centered matrix from PCA
		PCApipelineLite.setCorrectResultsForSTdevs(false);//not used, it uses centered matrix from PCA
		PCApipelineLite.setWriteFolder(pcaDir);
		PCApipelineLite.setLog2(false);//not used, it uses centered matrix from PCA
		PCApipelineLite.setDeSeqNorm(false);//not used, it uses centered matrix from PCA
		PCApipelineLite.setWriteFolderCorrected(correctedDir);
		
		//still have to fix this:
		String correctedFn = PCApipelineLite.getWriteFolderCorrected()+"PC_1-100.txt.gz";
		
//		LargeFileSorter fileSorter = new LargeFileSorter();
//		fileSorter.setFn(geneNameAdder.getWriteFnNoGeneSymbol());
//		fileSorter.setSortOrderFn(correctedFn);
//		fileSorter.setFnSorted(FileUtils.addBeforeExtention(geneNameAdder.getWriteFnNoGeneSymbol(), "sorted"));
//		
		FileSplitter fileSplitterCorrectedRatios = new FileSplitter();
		fileSplitterCorrectedRatios.setFn(correctedFn);
		fileSplitterCorrectedRatios.setGzipFiles(false);
		fileSplitterCorrectedRatios.setSplitRegex("__");
		fileSplitterCorrectedRatios.setSplitFirstColumnInstead(true);
		fileSplitterCorrectedRatios.setWriteFolder(correctedDir+"perGene/");
		
		// create exon expression files and calculate expression ratios per gene
		ExonExpressionMerger_InfiniteFileSizes exonExpressionMerger = new ExonExpressionMerger_InfiniteFileSizes();
		exonExpressionMerger.setJsonFN(this.getOutputFolder() + exonExpressionMerger.getNewJsonFN("SecondPass"));
		exonExpressionMerger.setSearchFolder(sTARv_SecondPass.getOutputRoot());
		exonExpressionMerger.setWriteFolder(sTARv_SecondPass.getOutputRoot());
		exonExpressionMerger.setWriteMergedFn(sTARv_SecondPass.getOutputRoot()+"expression_Exons.txt.gz");
		exonExpressionMerger.setEnsgToGeneSymbolFn(this.ensgToGeneSymbolFN);
		exonExpressionMerger.writeConfig();

		// add steps (in right order of course ;))
		List<Runnable> steps = new ArrayList<>();

		if (this.getPass() == 1 || this.getPass() == 2)
		{
			if (this.getFastQfilesFnFirstPass() == null)// if list of fastq files is not supplied yet, find them in the supplied folder(s)
				steps.add(this.fileSearcherFirstPass);
			steps.add(firstBuildVars);
			steps.add(sTARv_FirstPass);
			steps.add(laneMerger_StarOuptut_FirstPass);
			steps.add(mergeV_FirstPass);
		}
		if (this.getPass() == 2 || this.getPass() == -2)
		{
			if (this.getFastQfilesFnSecondPass() == null)// if list of fastq files is not supplied yet, find them in the supplied folder(s)
				steps.add(this.fileSearcherSecondPass);
			steps.add(secondBuildVars);
			steps.add(sTARv_SecondPass);
			steps.add(laneMerger_StarOuptut_SecondPass);
			steps.add(mergeV_SecondPass);
			steps.add(geneNameAdder);
			
			steps.add(fileSplitter);
			steps.add(ratioCalculator);
			steps.add(transpose);
			steps.add(centerPerRow);
			steps.add(correlationLarge);
			steps.add(pca);
			steps.add(averagesPerRow);
			steps.add(PCApipelineLite);
			//steps.add(fileSplitterCorrectedRatios);
			steps.add(exonExpressionMerger);
		}
		return steps;

	}

	/**Creates a new folder if it does not yet exist, unless the parent directory does not exist.
	 * 
	 * @param folderName name of the folder to be created.
	 */

	private void mkdir(String folderName)
	{
		if (folderName == null)
		{
			log("folderName is null: please set variable folderName\t");
			System.exit(2);
		}
		File folder = new File(folderName);
		File parent = new File(folder.getParent());
		if (!parent.exists())
		{
			log("Exiting, parent folder" + parent.getAbsolutePath() + " does not exist:\t");
			System.exit(2);
		}

		if (!new File(getOutputFolder()).exists())
			folder.mkdir();
	}

	private String getFastqFirstPass()
	{
		if (this.getFastQfilesFnFirstPass() != null)
			return this.getFastQfilesFnFirstPass();
		;
		return fileSearcherFirstPass.getWriteName();
	}

	private String getFastqSecondPass()
	{
		if (this.getFastQfilesFnSecondPass() != null)
			return this.getFastQfilesFnSecondPass();
		;
		return fileSearcherSecondPass.getWriteName();
	}

	public String includeWritePath(String name)
	{
		if (!name.contains("\\") && !name.contains("/"))
			return FileUtils.makeFolderNameEndWithSlash(this.outputFolder) + name;
		return name;
	}

	public static STAR_Pipeline readVars(	String jsonFilename,
											String executeClassName)// function should be in separate interface...
	{
		String jsonString = "";
		try
		{
			File file = new File(jsonFilename);
			if (!file.exists())
			{
				System.out.println("Json file does not exist:" + jsonFilename);
				System.exit(0);
			}
			jsonString = new String(Files.readAllBytes(Paths.get(jsonFilename)));
		} catch (Exception e)
		{
			e.printStackTrace();
		}

		Gson gson = new GsonBuilder().setPrettyPrinting().create();
		STAR_Pipeline var = gson.fromJson(	jsonString,
											STAR_Pipeline.class);
		var.jsonFN = jsonFilename;
		System.out.println(gson.toJson(var));
		return var;
	}

	public String getOutputFolder()
	{
		String outputFolder = FileUtils.makeFolderNameEndWithSlash(this.outputFolder);
		return outputFolder;
	}

	public STAR_Pipeline clone() throws CloneNotSupportedException
	{
		STAR_Pipeline copy = (STAR_Pipeline) super.clone();
		return copy;
	}

	public int getPass()
	{
		return pass;
	}

	private void setPass(int pass)
	{
		this.pass = pass;
	}

	public String getFastQfilesFnFirstPass()
	{
		return fastQfilesFirstPassFN;
	}

	public String getFastQfilesFnSecondPass()
	{
		return fastQfilesSecondPassFN;
	}

	private void setFastQfilesFN(String fastQfilesFN)
	{
		this.fastQfilesFirstPassFN = fastQfilesFN;
	}

	public void createConfig(String jsonFN)
	{
		this.slurmVars.jobType = new STAR_ClusterHandler();
		createConfig(	jsonFN,
						this,
						false,
						true);
	}

	@Override
	public HashMap<String, Integer> getStringAllowence()
	{//variables to be overwritten by child classes
		HashMap<String, Integer> stringAllowence = super.getStringAllowence();
		stringAllowence.put("writeNameComment",
							0);
		stringAllowence.put("writeName",
							0);
		stringAllowence.put("fastQFiles",
							0);
		stringAllowence.put("outputRoot",
							0);
		stringAllowence.putAll(this.slurmVars.getStringAllowence());
		stringAllowence.put("inputFolder_SpliceComment",
							0);
		stringAllowence.put("inputFolder_Splice",
							0);
		stringAllowence.put("writeFN_SpliceComment",
							0);
		stringAllowence.put("writeFN_Splice",
							0);
		stringAllowence.put("writeFolder_SplicePerGeneComment",
							0);
		stringAllowence.put("writeFolder_SplicePerGene",
							0);

		return stringAllowence;
	}

}
