package STAR;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.nio.file.Files;
import java.nio.file.Paths;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

public final class Variables implements Cloneable //A class that holds all the variables for the STAR package
{
	private int pass = 2;//1 for first pass only.		2 for first and second pass. 	-2 for secondpass only.
	private String jsonFN = null;//"E:/Groningen/Test/STAR/FastqMappingSTAR.json";
	private String outputFolder = null;//"E:/Groningen/Test/STAR/";
	private String inputFolder = null;//"J:/DATA/";
	private String fastQSearchStrings = null;//".fastq,.fq";
	private String forbiddenSearchStrings = null;//".md5";
	private String includeList = null;
	private String fastQfilesFN = null;//list of paths to fastQ files to include in the analysis (overwrites inputFolder, fastQSearchStrings and forbiddenSearchStrings)

	//Mapper variables
	private String outputRootSTAR = null;//outputFolder+"STAR/";
	private String genomeDir= null;//"/somedir/";
	private int batchSize = 20;
	//in case of star need to run:
	//STAR --runMode genomeGenerate --genomeDir /groups/umcg-gdio/tmp04/umcg-svandam/Data/RNAseq/Annotation/GCRh38/STARindex/  --genomeFastaFiles /groups/umcg-gdio/tmp04/umcg-svandam/Data/RNAseq/Annotation/GCRh38/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa --runThreadN 12
	private String version = null;//"STAR/2.5.1b-foss-2015b";
	private String threads = null;//"10";
	private String walltime = null;//"05:59:00";
	private String maxMemory = null;//"30gb";
	private String maxMemoryGenomeBuild = null;//"1gb";
	private String pairStrings = null;//"_R1_,_R2_,_1.fq,_2.fq";
	
	private String slurmUserName = null;//"umcg-svandam";
	private String finishedemailaddress = null;//"sipkovandam@gmail.com";
	private String sjdbGTFfile = null;
	private String STAR_Arguments = null;//--clip5pNbases 102 --outFilterMismatchNmax "999" --outFilterMismatchNoverLmax 0.04
	private String keepBAMsContaining = "RNA14-00231_S1,RNA14-00254_S7,RNA14-00258_S4";//if a bam file contains any of these strings it is kept
	private transient String readFN_Splice = null;//in the second pass this file is used to retrieve the splice sites from (if the first pass is also ran, this becomes the output file of the first pass)
	
	//FeatureCounts
	private String featureCounts = null;
	private String featureCountsOptions = null;//"-T 25 --largestOverlap -f -t exon -O"-f -t exon options make sure it counts and reports reads per exon instead of gene ,-O (capital o) option makes sure reads are counted whenever it overlaps multiple genes
	
	//SpliceSites variables
	private String inputFolder_Splice = outputRootSTAR;
	private transient String writeFN_Splice = null;//file where all the splice variants merged into 1 file with info on number of spliced reads overlapping each in all samples together
	private String annotationFN = null;//"E:/Groningen/Data/Annotation/GRCh38/ChromosomeInfo/GenePositionInfo_withMT.txt";
	private double minPercentage = 0.1;//number of samples this splice variant has to be present in for it to be included in every file (thus also in those files where it is not expressed)
	private int readCutoff = 8;
	private String excludeFromReferenceSamplesContaining = null;//"160613_SN163_0713_AC8NKTACXX_L5_CAACTA";//cancer sample (BRCA2) //"0667_AC7CNEACXX,0694_AC8N1LACXX,0713_AC8NKTACXX";//all clinical samples including parents  
	private double minPercentageSamplesPresentForSpliceToBeWrittenToOwnFile = 0.25;//Minimun percentage of files in which a splice variant has to be expressed in order to be written to its own file
	private String writeFolder_Splice = null;//if null becomes: new File(var.writeFN).getParent();

	//SpliceSitesPerGene variables
	private String ensgToGeneSymbolFN = null;//"E:/Groningen/Data/Annotation/GRCh38/EnsgToGeneSymbol.txt";//"/groups/umcg-gdio/tmp04/umcg-svandam/Data/RNAseq/Annotation/GCRh38/EnsgToGeneSymbol.txt";//
	
	public Variables()
	{
		
	}
	
	
	public void writeVars()
	{
		if(this.jsonFN == null)
			this.jsonFN = new File(getWriteFN_Splice()).getParent()+"/SpliceSites.json";
		writeVars(this.jsonFN);
		
	}
	void writeVars(String jsonFN)
	{
		Gson gson = new GsonBuilder().serializeNulls().setPrettyPrinting().create();
		//System.out.println(gson.toJson(this));
		//System.out.println(gson.toString());
		System.out.println(getWritePath(jsonFN));
		this.jsonFN = jsonFN;
		try
		{
			BufferedWriter writer = new BufferedWriter(new FileWriter(new File(getWritePath(jsonFN))));
			
			writer.write(gson.toJson(this));
			writer.close();
		}catch(Exception e){e.printStackTrace();}
	}

	public String getFolderName(String fn) 
	{
		if(!fn.endsWith("/") && !fn.endsWith("\\"))
			fn = fn+"/";
		return fn;
	}

	public String getWritePath(String name)
	{
		if(!name.contains("\\") && !name.contains("/"))
			return getFolderName(this.outputFolder)+name;
		return name;
	}
	static Variables readVars(String jsonFilename)
	{
		String jsonString = "";
		try
		{
			File file = new File(jsonFilename);
			if(!file.exists())
			{
				System.out.println("Json file does not exist:" + jsonFilename);
				System.exit(0);
			}
			jsonString = new String(Files.readAllBytes(Paths.get(jsonFilename)));
		}catch(Exception e){e.printStackTrace();}
		
		Gson gson = new GsonBuilder().setPrettyPrinting().create();
		Variables var = gson.fromJson(jsonString, Variables.class);

		var.setWriteFolder_Splice();

		var.jsonFN=jsonFilename;
		System.out.println(gson.toJson(var));
		return var;
	}

	public String getWriteFolder_Splice() {
		setWriteFolder_Splice();
		return writeFolder_Splice;
	}
	public String getJsonFN() {
		return jsonFN;
	}


	public String getOutputFolder() {
		return outputFolder;
	}


	public String getFastQSearchStrings() {
		return fastQSearchStrings;
	}


	public String getForbiddenSearchStrings() {
		return forbiddenSearchStrings;
	}


	public String getOutputRootSTAR() {
		return outputRootSTAR;
	}
	public void setOutputRootSTAR(String newOutputRoot) {
		outputRootSTAR = newOutputRoot;
	}

	public String getGenomeDir() {
		return genomeDir;
	}
	public void setGenomeDir(String newDir) {
		genomeDir = newDir;
	}

	public String getVersion() {
		return version;
	}


	public String getThreads() {
		return threads;
	}


	public String getWalltime() {
		return walltime;
	}


	public String getMaxMemory() {
		return maxMemory;
	}


	public String getPairStrings() {
		return pairStrings;
	}


	public String getSlurmUserName() {
		return slurmUserName;
	}


	public String getFinishedemailaddress() {
		return finishedemailaddress;
	}


	public String getSjdbGTFfile() {
		return sjdbGTFfile;
	}


	public String getInputFolder_Splice() {
		if(inputFolder_Splice == null)
			inputFolder_Splice=outputRootSTAR;
		return inputFolder_Splice;
	}


	public String getWriteFN_Splice() {
		if(writeFN_Splice==null)
			writeFN_Splice = outputRootSTAR+"SJ_Merged.out.tab";
		return writeFN_Splice;
	}
	
	public void setWriteFN_Splice(String newWriteFN_Splice) {
		writeFN_Splice=newWriteFN_Splice;
	}


	public String getAnnotationFN() {
		return annotationFN;
	}

	public double getMinPercentage() {
		return minPercentage;
	}

	public String getEnsgToGeneSymbolFN() {
		return ensgToGeneSymbolFN;
	}

	private void setWriteFolder_Splice() {
		if(this.writeFolder_Splice ==null)
			this.writeFolder_Splice = new File(getWriteFN_Splice()).getParent()+"/";
	}

	public String getInputFolder() {
		return inputFolder;
	}
	
	public Variables clone() throws CloneNotSupportedException {
	      return (Variables) super.clone();
	}


	public void setInputFolder_Splice(String inputFolder) {
		this.inputFolder_Splice = inputFolder;
	}


	/**
	 * @return the readCutoff
	 */
	public int getReadCutoff() {
		return readCutoff;
	}


	/**
	 * @param readCutoff the readCutoff to set
	 */
	public void setReadCutoff(int readCutoff) {
		this.readCutoff = readCutoff;
	}


	/**
	 * @return the excludeSamplesCalculatingStatsContaining
	 */
	public String getExcludeFromReferenceSamplesContaining() {
		return excludeFromReferenceSamplesContaining;
	}


	/**
	 * @param excludeSamplesCalculatingStatsContaining the excludeSamplesCalculatingStatsContaining to set
	 */
	public void setExcludeFromReferenceSamplesContaining(String excludeSamplesCalculatingStatsContaining) {
		this.excludeFromReferenceSamplesContaining = excludeSamplesCalculatingStatsContaining;
	}


	/**
	 * @return the minSamplesForSpliceToBeWrittenToOwnFile
	 */
	public double getMinPercentageSamplesPresentForSpliceToBeWrittenToOwnFile() {
		return minPercentageSamplesPresentForSpliceToBeWrittenToOwnFile;
	}


	public String getSTAR_Arguments() {
		return STAR_Arguments;
	}


	public void setSTAR_Arguments(String sTAR_Arguments) {
		STAR_Arguments = sTAR_Arguments;
	}


	public int getBatchSize() {
		return batchSize;
	}


	public void setBatchSize(int batchSize) {
		this.batchSize = batchSize;
	}


	public int getPass() {
		return pass;
	}


	public String getFeatureCounts() {
		return featureCounts;
	}


	public String getFeatureCountsOptions() {
		return featureCountsOptions;
	}


	public void setFeatureCountsOptions(String featureCountsOptions) {
		this.featureCountsOptions = featureCountsOptions;
	}


	public String getKeepBAMsContaining() {
		return keepBAMsContaining;
	}


	public void setKeepBAMsContaining(String keepBAMsContaining) {
		this.keepBAMsContaining = keepBAMsContaining;
	}


	public String getReadFN_Splice() {
		return readFN_Splice;
	}


	public void setReadFN_Splice(String readFN_Splice) {
		this.readFN_Splice = readFN_Splice;
	}


	public String getIncludeList() {
		return includeList;
	}


	public void setIncludeList(String includeList) {
		this.includeList = includeList;
	}


	public String getFastQfilesFN() {
		return fastQfilesFN;
	}


	public String getMaxMemoryGenomeBuild() {
		return maxMemoryGenomeBuild;
	}


	public void setMaxMemoryGenomeBuild(String maxMemoryGenomeBuild) {
		this.maxMemoryGenomeBuild = maxMemoryGenomeBuild;
	}

//	public boolean filePathsExist()
//	{
//		File[] files = new File[7]; 
//		
//		for(File file : files)
//		{
//			if(file.exists())
//				continue;
//			System.out.println("THIS FILE/FOLDER DOES NOT EXIST!\n" + file.getAbsolutePath());
//			return false;
//		}
//		return true;
//	}
	
}
