package Kallisto;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;

import org.junit.Rule;
import org.junit.rules.TemporaryFolder;

import Kallisto.FastQtoExpression.Vars;
import Tools.ExecCommand;
import Tools.FileUtils;
import Tools.JSONutil;
import Tools.Util;

public class KallistoSlurm {

	@Rule
	public static TemporaryFolder tempFolder = new TemporaryFolder();

	static String fastQFiles = "E:/Groningen/Test/JSON/fastQfiles.txt";
	static String writeFolder = null;
	static String scriptsFolderName = null;
	static String logsFolder = null;
	static String errorsFolder = null;
	static String kallistoVersion = "Kallisto/0.42.2.1-goolf-1.7.20";//"Kallisto/0.42.4-goolf-1.7.20";
	static String kallistoIndexFile = "/groups/umcg-wijmenga/tmp04/umcg-svandam/Data/RNAseq/Annotation/hg19.v75.cdna.all.42.2.idx";
	static int kallistoThreads = 4;
	static String kallistoWalltime = "05:59:00";
	static String kallistoMaxMemory = "8gb";
	static String kallistoOutputRoot = null;
	static String mappingPercentagesFN = null;
	
	static String slurmUserName = "umcg-svandam";
	static String finishedemailaddress = "sipkovandam@gmail.com";
	static String tsvFilesToShScriptFN = null;
	
	// pairStrings to be in pairs. String in uneven index is replaced with that in the even index to get the paired end partner.
	//All file names containing a string in one of the uneven indexes is skipped)
	static String[] pairStrings = new String[]{"_R1_","_R2_","_1.fq","_2.fq"};
	
	public static void main(String[] args) throws Exception
	{
		checkArgs(args);
		if(writeFolder == null)
			writeFolder = new File(fastQFiles).getParent()+"/";
		System.out.println(writeFolder);
		if(scriptsFolderName == null)
			scriptsFolderName = writeFolder+"Scripts/";
		if(logsFolder == null)
			logsFolder = writeFolder+"Logs/";
		if(errorsFolder == null)
			errorsFolder = writeFolder+"Errors/";
		if(kallistoOutputRoot == null)
			kallistoOutputRoot = writeFolder+"Results/";
		if(mappingPercentagesFN == null)
			mappingPercentagesFN = new File(kallistoOutputRoot).getParent()+"/mappingPerSample.txt";
		if(tsvFilesToShScriptFN == null)
			tsvFilesToShScriptFN = new File(scriptsFolderName).getParent()+"/scriptNumberToFiles.txt";
		System.out.println("scriptsFolderName=" + scriptsFolderName);
		
		ArrayList<String> fastQ_Files = FileUtils.readArrayList(fastQFiles);  
		FileUtils.makeDir(scriptsFolderName);
		FileUtils.makeDir(logsFolder);
		FileUtils.makeDir(errorsFolder);
		FileUtils.makeDir(kallistoOutputRoot);
		
		createSlurmFiles(fastQ_Files);
		runKallistoSlurmScripts();
	}

	private static void runKallistoSlurmScripts() throws FileNotFoundException, IOException 
	{
		//carefull with this. Kallisto200sh script does not update! You need to delete and reimport project of that!
		tempFolder.create();
		String kallistoShellFN = "Kallisto200.sh";

		System.out.println("exists = " + new File(kallistoShellFN).exists());
		String command = "bash " + kallistoShellFN + " " + scriptsFolderName +"*.sh "+kallistoOutputRoot+" "+ slurmUserName + " "+ mappingPercentagesFN +" " + finishedemailaddress;
		System.out.println("Shellcommand = " + command);
		ExecCommand exec = new ExecCommand(command);
		//System.out.println("execute output: \n"+exec.getOutput());
		System.out.println("execute error: \n"+exec.getError());
	}
//  Redundant, Jesse heeft hier voor een paar classes gegeven
//	private static void runShell(String command) {
//		try
//        {            
//            Runtime rt = Runtime.getRuntime();
//            Process proc = rt.exec(command);
//            InputStream stderr = proc.getErrorStream();
//            InputStreamReader isr = new InputStreamReader(stderr);
//            BufferedReader br = new BufferedReader(isr);
//            String line = null;
//            while ( (line = br.readLine()) != null)
//                System.out.println(line);
//            int exitVal = proc.waitFor();
//            System.out.println("Process exitValue: " + exitVal);
//        } catch (Throwable t)
//        {
//        	t.printStackTrace();
//        }
//	}

	private static void createSlurmFiles(ArrayList<String> fastQ_Files) throws Exception {
		int fileNumber =1;
		BufferedWriter tsvFilenameWriter = FileUtils.createWriter(tsvFilesToShScriptFN);
		out: for(int f = 0; f < fastQ_Files.size(); f++)
		{
			String shellFN = scriptsFolderName+fileNumber+".sh";
			
			String fastqFN = fastQ_Files.get(f);
			
			//continue if it is the second file of a paired end sequenced sample
			String pairedStringForward = null;
			String pairedStringReverse = null;
			for(int p = 1; p < pairStrings.length; p+=2)
				if(fastqFN.contains(pairStrings[p]))
					continue out;
			
			for(int p = 0; p < pairStrings.length; p+=2)
				if(fastqFN.contains(pairStrings[p]))
				{
					pairedStringForward = pairStrings[p];
					pairedStringReverse = pairStrings[p+1];
				}
			
			BufferedWriter writer = FileUtils.createWriter(shellFN);
			writeSlurmCommands(writer, f, fastqFN, pairedStringForward, pairedStringReverse,fileNumber,tsvFilenameWriter);
			System.out.println("File written to: " + shellFN);
			writer.close();
			fileNumber++;
		}
		tsvFilenameWriter.close();
	}

	private static void writeSlurmCommands(BufferedWriter writer, int lineNum, String fn, String pairedStringForward, String pairedStringReverse, int fileNumber, BufferedWriter tsvFilenameWriter) throws Exception {
		// TODO Auto-generated method stub
		writer.write("#!/bin/bash"+"\n");
		writer.write("#SBATCH --job-name=kallisto"+fileNumber+"\n");
		writer.write("#SBATCH --output="+logsFolder+fileNumber+".out\n");
		writer.write("#SBATCH --error="+errorsFolder+fileNumber+".err\n");
//		writer.write("#SBATCH --partition=medium\n");
		writer.write("#SBATCH --time="+kallistoWalltime+"\n");
		writer.write("#SBATCH --cpus-per-task "+kallistoThreads+"\n");
		writer.write("#SBATCH --mem "+kallistoMaxMemory+"\n");
		writer.write("#SBATCH --nodes 1\n");
		writer.write("#SBATCH --export=NONE\n");
		writer.write("#SBATCH --get-user-env=L\n");
//		writer.write("#SBATCH --open-mode=append\n");
		writer.write("module load "+kallistoVersion+"\n");
		writer.write("kallisto version\n");

		writeKallisto(writer, fn, pairedStringForward, pairedStringReverse, fileNumber, tsvFilenameWriter);
	}
	private static void writeKallisto(BufferedWriter writer, String fn, String pairedStringForward, String pairedStringReverse, int fileNumber, BufferedWriter tsvFilenameWriter) throws Exception {

		File file = new File(fn);
		boolean iris = false;
		if(fn.contains("Iris"))//does not do anything atm (Iris changed the structure of her folders)
			iris = true;
		if(pairedStringForward != null && file.getName().contains(pairedStringForward))
			writePairedEndCommand(file, pairedStringForward, pairedStringReverse, iris, writer, fileNumber, tsvFilenameWriter);
		else if(file.getName().contains(".fq")) //single end
			writeSingleEndCommand(file, iris, writer, fileNumber, tsvFilenameWriter);
		writer.close();
	}

	private static void writePairedEndCommand(File file, String pairedStringForward, CharSequence pairedStringReverse, boolean iris, BufferedWriter writer, int fileNumber, BufferedWriter tsvFilenameWriter) throws IOException {
		String outputFolder = file.getName().split(pairedStringForward)[0]+"/";
//		if(iris)//Iris samples need a special treatment because the name of the file is 1 folder higher in the folder hierarchie ><
//		{
//			String[] parentFolders = file.getParent().split("\\\\");
//			outputFolder = file.getName().replace("_1.fq.gz", "")+"/";
//		}
		writer.write("mkdir " + kallistoOutputRoot + outputFolder+"\n");
		String fileName = "\""+file.getName()+"\"";
		fileName = fileName.replace("\\", "/");
		String fastqWithPath = file.getPath().replace("\\", "/");
		String line = 	"kallisto quant --bias -t "+kallistoThreads+" -i " + kallistoIndexFile + " -o " + kallistoOutputRoot + outputFolder + " \"" + 
						fastqWithPath + "\" \"" + fastqWithPath.replace(pairedStringForward, pairedStringReverse) + "\" &> "+
						kallistoOutputRoot + outputFolder +outputFolder.replace("/","")+".err" +"\n";
		writer.write(line);
		tsvFilenameWriter.write(kallistoOutputRoot+outputFolder+"abundance.tsv"+"\t"+fileNumber+".sh"+"\t"+fastqWithPath+"\n");
	}

	private static void writeSingleEndCommand(File file, boolean iris, BufferedWriter writer, int fileNumber, BufferedWriter tsvFilenameWriter) throws IOException 
	{
		String outputFolder = file.getName().replace(".fq.gz","/");
//		if(iris)//Iris samples need a special treatment because the name of the file is 1 folder higher in the folder hierarchie ><
//		{
//			String[] parentFolders = file.getParent().split("\\\\");
//			outputFolder = parentFolders[parentFolders.length-1]+"/";
//		}
		writer.write("mkdir " + kallistoOutputRoot + outputFolder+"\n");
		String fileName = "\""+file.getName()+"\"";
		fileName = fileName.replace("\\", "/");
		String fastqWithPath = file.getPath().replace("\\", "/");
		String line = 	"kallisto quant --bias -t "+kallistoThreads+" --single --fragment-length=200 --sd=20 -i " + kallistoIndexFile + " -o " + kallistoOutputRoot + 
						outputFolder + " \"" + fastqWithPath + "\" &> "+ kallistoOutputRoot + outputFolder+ 
						"kallisto_"+outputFolder.replace("/","")+".err" +"\n";
		writer.write(line);
		tsvFilenameWriter.write(kallistoOutputRoot+outputFolder+"abundance.tsv"+"\t"+fileNumber+".sh"+"\t"+fastqWithPath+"\n");
	}
	static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko") && args.length < 1)
			return;
		if(args.length < 1)
		{
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

		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase())
			{
//				case "fastQFiles":
//					var.JSON_FN =value;
//					var = new JSONutil<Vars>().read(var.JSON_FN, var);
//					break;
				case "fastqfiles":
					fastQFiles = value;
					break;
				case "kallistoindexfile":
					kallistoIndexFile = value;
					break;
				case "scriptsfoldername":
					scriptsFolderName = value;
					break;
				case "logsfolder":
					logsFolder = value;
					break;
				case "errorsfolder":
					errorsFolder = value;
					break;
				case "kallistooutputroot":
					kallistoOutputRoot = value;
					break;
				case "kallistoversion":
					kallistoVersion = value;
					break;
				case "kallistothreads":
					kallistoThreads = Integer.parseInt(value);
					break;
				case "kallistowalltime":
					kallistoWalltime = value;
					break;
				case "kallistomaxmemory":
					kallistoMaxMemory = value;
					break;
				case "pairstrings":
					pairStrings = value.split(",");
					break;
				case "slurmusername":
					slurmUserName = value;
					break;
				case "finishedemailaddress":
					finishedemailaddress = value;
					break;
				case "writefolder":
					writeFolder = value;
					break;
				case "mappingpercentagesfn":
					mappingPercentagesFN = value;
					break;
				case "tsvfilestoshscriptfn":
					tsvFilesToShScriptFN = value;
					break;
				default:
					System.out.println("Incorrect argument supplied:\n"+ args[a] +"\nexiting");
					System.exit(1);
			}
		}
	}
}
