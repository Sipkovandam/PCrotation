package STAR;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Hashtable;
import java.util.concurrent.atomic.AtomicInteger;

import com.sun.msv.datatype.xsd.datetime.Util;

import Kallisto.Slurm;
import Tools.FileUtils;
import Tools.JSONutil;
import Tools.SearchFilesInDirectories;

public class FastqMappingSTAR 
{
	//This class does the following:
	//1. Finds all fastQ files in a folder and subfolders
	//2. Determines which fastQ files are pairs
	//3. Maps the files using Kallisto
	//4. Creates a file describing the number of reads per sample mapped
	//5. Sums reads of transcritps to genes
	//6. Creates 1 large matrix for genes and one for transcripts
	//7. Creates 1 large matrix for samples >70% mapping
	
	//groups/umcg-wijmenga/tmp04/umcg-svandam/FeatureCounts/subread-1.5.1-source/bin/featureCounts
	
	public static void main(String args[]) throws Exception
	{
		//args = new String[]{"json=E:/Groningen/Test/STAR/var.json"};
		//args = new String[]{"json=E:/Groningen/Splicing/100BPcap_analysis/STAR_config_200Samples.json"};
		Variables v = checkArgs(args);
		v.writeVars();
		
		if(!new File(v.getOutputFolder()).exists())
			new File(v.getOutputFolder()).mkdir();
		
		final String fastQFiles=findFastQFiles(v);
		
		Variables v_FirstPassMerge = v.clone();
		v_FirstPassMerge.setWriteFN_Splice(v_FirstPassMerge.getOutputFolder()+"SJ_Merged_1stPass.out.tab");
		if(v.getPass() == 1 || v.getPass() == 2)
		{
			//first pass
			System.out.println("Running first pass STAR");
			v_FirstPassMerge.setOutputRootSTAR(v_FirstPassMerge.getOutputFolder()+"Results_1stPass/");
			v_FirstPassMerge.setInputFolder_Splice(new File(v_FirstPassMerge.getOutputRootSTAR())+"/Results/");
			//STAR_Slurm(v_FirstPassMerge, fastQFiles, "None","false");		
			System.out.println("merging splice files from 1st pass");
			
			mergeSpliceFiles(v_FirstPassMerge);
		}
		
		if(v.getPass() == 2 || v.getPass()==-2)
		{
			//rebuild the genome.
			Variables buildVars = v.clone();
			buildVars.setOutputRootSTAR(buildVars.getOutputFolder()+"reBuild/");
			buildVars.setSTAR_Arguments(buildVars.getSTAR_Arguments()+" --sjdbInsertSave All");
			buildVars.setReadFN_Splice(v_FirstPassMerge.getWriteFN_Splice());
			String genomeDirWithSplice=generateNewGenome(buildVars,fastQFiles);
			
			System.out.println("genomeDirWithSplice=" + genomeDirWithSplice);
			
			//second pass
			System.out.println("Running 2nd pass STAR");
			Variables v_SecondPassMerge = v.clone();
			v_SecondPassMerge.setGenomeDir(genomeDirWithSplice);
			v_SecondPassMerge.setOutputRootSTAR(v_SecondPassMerge.getOutputFolder()+"Results_2ndPass/");
			v_SecondPassMerge.setInputFolder_Splice(new File(v_SecondPassMerge.getOutputRootSTAR())+"Results/");
			STAR_Slurm(v_SecondPassMerge, fastQFiles, "BAM Unsorted","true");	//last argument indicates whether STAR should also determine expression per gene (this makes it quite a lot slower I think)
			v_SecondPassMerge.setWriteFN_Splice(v_SecondPassMerge.getOutputFolder()+"SJ_Merged_2ndpass.out.tab");
			mergeSpliceFiles(v_SecondPassMerge);
			System.out.println("Second pass output files written to:" + v_SecondPassMerge.getOutputRootSTAR());
		}
	}
	private static String generateNewGenome(Variables v, String fastQFiles) throws Exception {
		String genomeDirWithSplice=v.getOutputRootSTAR()+"Results/_STARgenome/";//STAR automatically adds "_STARgenome/"
		STAR_Slurm(v, fastQFiles, "None","true", v.getOutputRootSTAR());
		return genomeDirWithSplice;
	}
	private static void mergeSpliceFiles(Variables v_FirstPassMerge) throws Exception {
		SpliceSites.run(v_FirstPassMerge);
	}
	
	private static void addToHash(String file, Hashtable<String, int[]> spliceHash){
		BufferedReader readerSpliceFile;
		try {
			readerSpliceFile = FileUtils.createReader(file);
		
			readerSpliceFile.lines().forEach(line -> addToSpliceHash(line,spliceHash));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private static void addToSpliceHash(String line, Hashtable<String, int[]> spliceHash) 
	{
		//get the bit that defines the splice variant
		AtomicInteger tab = new AtomicInteger(0);
		AtomicInteger pos = new AtomicInteger(0);
		AtomicInteger splitChar = new AtomicInteger(0);
		line.chars().forEach(c -> {
			if(c== 9)//9 == tab
				tab.getAndIncrement();
				if(tab.get()==6)
					splitChar.set(pos.get());
				pos.getAndIncrement();
					});
		String spliceVar = line.substring(0,splitChar.get());
		
		String[] eles = line.split("\t");

		int observations = Integer.parseInt(eles[7]);
		int maxOverhang = Integer.parseInt(eles[9]);
		
		if(!spliceHash.containsKey(spliceVar)) 
		{
			int [] values = new int[4];
			values[0]=observations;//how many reads overlap this junction
			values[1]=0;//how often it is observed in multimapped reads (multimapped reads are ignored using the settings I use)
			values[2]=maxOverhang;//max overhang
			values[3]=1;//number of samples it is observed in
			spliceHash.put(spliceVar, values);
		}
		else
		{
			int[] values = spliceHash.get(line);
			values[0]+=observations;//how many reads overlap this junction
			values[1]+=0;//how often it is observed in multimapped reads (multimapped reads are ignored using the settings I use)
			if(maxOverhang>values[2])
				values[2]=maxOverhang;//max overhang
			values[3]+=1;//number of samples it is observed in
			spliceHash.put(line, values);
		}
	}
	private static void STAR_Slurm(Variables v, String fastQFiles) throws Exception
	{
		STAR_Slurm(v, fastQFiles, "Full", "false");
	}
	private static void STAR_Slurm(Variables v, String fastQFiles, String outputFormat, String countExpression) throws Exception 
	{
		STAR_Slurm(v, fastQFiles, outputFormat, countExpression,  null); 
	}

	private static void STAR_Slurm(Variables v, String fastQFiles, String outputFormat, String countExpression, String buildDir) throws Exception {
		new Slurm().run(new String[]{	"fastQFiles="+fastQFiles,
										"mapper="+v.getVersion(),
										"writefolder="+v.getOutputRootSTAR(),
										"threads="+v.getThreads(),
										"walltime="+v.getWalltime(),
										"maxMemory="+v.getMaxMemory(),
										"pairStrings="+v.getPairStrings(),
										"slurmUserName="+v.getSlurmUserName(),
										"genomeDir="+v.getGenomeDir(),
										"GTFfile="+v.getSjdbGTFfile(),
										"sjdbFileChrStartEnd="+v.getReadFN_Splice(),
										"finishedemailaddress="+v.getFinishedemailaddress(),
										"outMode="+outputFormat,
										"countExpression="+countExpression,
										"STAR_Arguments="+v.getSTAR_Arguments(),
										"batchsize="+v.getBatchSize(),
										"buildDir="+buildDir,
										"featurecounts="+v.getFeatureCounts(),
										"featurecountsoptions="+v.getFeatureCountsOptions(),
										"keepBAMsContaining="+v.getKeepBAMsContaining(),
										"maxMemoryGenomeBuild="+v.getMaxMemoryGenomeBuild()});
		
	}
	
	private static String findFastQFiles(Variables v) throws IOException 
	{ 
		String fastQFiles = null;
		if(v.getFastQfilesFN() != null) 
			return v.getFastQfilesFN();
		
		fastQFiles = v.getOutputFolder()+"fastQfiles.txt";
		SearchFilesInDirectories.main(new String[]{	"foldername="+v.getInputFolder(),
													"writefn="+fastQFiles,
													"searchstrings="+v.getFastQSearchStrings(),
													"forbiddenstrings=" + v.getForbiddenSearchStrings(),
													"requiredstringfn=" + v.getIncludeList()});
		return fastQFiles;
	}
	
	static Variables checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko") && args.length < 1)
			return new Variables();
		if(args.length < 1)
		{
			System.out.println("Script requires the following argumetns:\n"
					+ "1. filename=<expressionFN.txt> - Expression file with genes on rows samples on columns\n"
					+ "2. writeFolder=<writeFolderFN.txt> - Folder where the files will be written (default=parentFolder(input.txt))\n"
					+ "3. geoFN=<geoFn.txt> - Optional file with geometric mean per gene to use (default=null)\n");
			System.exit(1);
		}

		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase())
			{
				case "json":
					return Variables.readVars(value);
				default:
					System.out.println("Incorrect argument supplied:\n"+ args[a] +"\nexiting");
					System.exit(1);
			}
		}
		return null;
	}
}
