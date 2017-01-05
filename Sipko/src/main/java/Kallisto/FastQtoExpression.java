package Kallisto;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import Tools.FileUtils;
import Tools.JSONutil;
import Tools.SearchFilesInDirectories;

public class FastQtoExpression 
{
	//This class does the following:
	//1. Finds all fastQ files in a folder and subfolders
	//2. Determines which fastQ files are pairs
	//3. Maps the files using Kallisto
	//4. Creates a file describing the number of reads per sample mapped
	//5. Sums reads of transcritps to genes
	//6. Creates 1 large matrix for genes and one for transcripts
	//7. Creates 1 large matrix for samples >70% mapping
	
	public static class Vars
	{
		String outputFolder = "E:/Groningen/Test/JSON/";
		String JSON_FN = "E:/Groningen/Test/JSON/var.json";
		String inputFolder = "J:/DATA/";
		String fastQSearchStrings = ".fastq,.fq";
		String forbiddenSearchStrings = ".md5";
		double mappingCutoff = 0.7;
		
		//Kallisto variables
		String kallistoOutputRoot = outputFolder+"Kallisto/";
		String kallistoVersion = "Kallisto/0.42.2.1-goolf-1.7.20";//"Kallisto/0.42.4-goolf-1.7.20";
		String kallistoIndexFile = "/groups/umcg-wijmenga/tmp04/umcg-svandam/Data/RNAseq/Annotation/hg19.v75.cdna.all.42.2.idx";
		String kallistoThreads = "4";
		String kallistoWalltime = "05:59:00";
		String kallistoMaxMemory = "8gb";
		String pairStrings = "_R1_,_R2_,_1.fq,_2.fq";
		
		String slurmUserName = "umcg-svandam";
		String finishedemailaddress = "sipkovandam@gmail.com";
		
		//file merging Kallisto Files, creating expression per gene and creating cutoff files variables
		String kallistoOutputColumn = "2";//2 for counts, 3 for tpm values
		String kallistoThreshold = "0.7";
		String transcriptsToGenesFN = "";
		String minpercentagefeaturesexpressed = "10";

	}
	static Vars var = new Vars();
	
	public static void main(String args[]) throws Exception
	{
		//args = new String[]{"json=E:/Groningen/Test/JSON/var.json"};
		checkArgs(args);
		new JSONutil<Vars>().write(var.JSON_FN, var);
		
		String fastQFiles = findFastQFiles();
		
		kallistoSlurm(fastQFiles);		
		
		combineKallistoOutput();
	}
	
	private static void combineKallistoOutput() throws Exception {
		String mappingPercentagesFN = var.kallistoOutputRoot+"mappingPerSample.txt";
		String tsvFilesToShScriptFN = var.kallistoOutputRoot+"scriptNumberToFiles.txt";
		CombineKallisto.main(new String[]{	"kallistoOutputFolder="+var.kallistoOutputRoot,
											"kallistoColumn="+var.kallistoOutputColumn,
											"transcriptsToGenesFN="+var.transcriptsToGenesFN,
											"threshold="+var.kallistoThreshold,
											"minpercentagefeaturesexpressed="+var.minpercentagefeaturesexpressed,
											"mappingpercentagesfn="+mappingPercentagesFN,
											"tsvthatshouldbetherefn="+tsvFilesToShScriptFN});
	}

	private static void kallistoSlurm(String fastQFiles) throws Exception {
		new Slurm().run(new String[]{"fastQFiles="+fastQFiles,
										"kallistoIndexFile="+var.kallistoIndexFile,
										"writefolder="+var.kallistoOutputRoot,
										"kallistoVersion="+var.kallistoVersion,
										"kallistoThreads="+var.kallistoThreads,
										"kallistoWalltime="+var.kallistoWalltime,
										"kallistoMaxMemory="+var.kallistoMaxMemory,
										"pairStrings="+var.pairStrings,
										"slurmUserName="+var.slurmUserName,
										"finishedemailaddress="+var.finishedemailaddress});
		
	}
	private static String findFastQFiles() throws IOException 
	{ 
		String fastQFiles = var.outputFolder+"fastQfiles.txt";
		SearchFilesInDirectories.main(new String[]{	"foldername="+var.inputFolder,
													"writefn="+fastQFiles,
													"searchstrings="+var.fastQSearchStrings,
													"forbiddenstrings=" + var.forbiddenSearchStrings});
		return fastQFiles;
	}
	static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko") && args.length < 1)
			return;
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
					var.JSON_FN =value;
					var = new JSONutil<Vars>().read(var.JSON_FN, var);
					break;
				default:
					System.out.println("Incorrect argument supplied:\n"+ args[a] +"\nexiting");
					System.exit(1);
			}
		}
	}
}
