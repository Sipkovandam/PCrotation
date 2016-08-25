package PrepData;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class RunKallistoOnFastqFiles 
{
	public static void main(String[] args) throws Exception
	{
//		String folderFN = "E:/Groningen/Data/Monogenetic_disease_samples_RadBoud/";
//		String rootName = "/Volumes/Promise_RAID/GeneNetwork/Sipko/Monogenetic_disease_samples_RadBoud/";
//		String indexFile = "/Volumes/Promise_RAID/juha/PublicRNASeq/kallisto/data/hg19.v75.cdna.all.42.2.idx";
//		String scriptName = "E:/Groningen/Data/Monogenetic_disease_samples_RadBoud/KallistoRunControls.sh";
		
		String folderFN = "E:/Groningen/Data/PublicSamples/07-2016/5GPM_samples/5gpmFastQs.txt";
		//output folder on the server
		String rootName = "/groups/umcg-gcc/tmp04/umcg-rkanninga/160613_SN163_0713_AC8NKTACXX/";
		String indexFile = "/target/gpfs2/groups/umcg-wijmenga/tmp02/projects/Celiac_Disease/gene_expression/RNASEQ/RNAseq_NG_data_08032016/HiSeq/stimulated_gluten_specific_Tcell_clones_TCC23052016/Kallisto/hg19.v75.cdna.all.42.2.idx";
		String scriptName = "E:/Groningen/Data/PublicSamples/07-2016/5GPM_samples/KallistoRun5gpmAll.sh";
		boolean iris = true; //takes the folder name for sample/column names in the files
		
		String PairedStringForward = "_1.fq";
		String PairedStringReverse = "_2.fq";
		
		createShellScript(folderFN, scriptName, PairedStringReverse, PairedStringForward, iris, rootName, indexFile);
	}
	private static void createShellScript(String folderFN, String scriptName, String PairedStringReverse, String PairedStringForward, String rootName, String indexFile) throws Exception 
	{
		createShellScript(folderFN, scriptName, PairedStringReverse, PairedStringForward, rootName, indexFile);
	}
	private static void createShellScript(String folderFN, String scriptName, String PairedStringReverse, String PairedStringForward, boolean iris, String rootName, String indexFile) throws Exception {
		File[] files = null;
		if(folderFN.contains(".txt"))
			files = getFileNames(folderFN);
		else
			files = new File(folderFN).listFiles();
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(scriptName)));
		writer.write("ml Kallisto\n");
		for(File f : files)
		{
			if(f.getName().contains(PairedStringReverse))
				continue;
			else if(f.getName().contains(PairedStringForward))
			{
				String outputFolder = f.getName().split(PairedStringForward)[0]+"/";
				if(iris)
				{
					String[] parentFolders = f.getParent().split("\\\\");
					outputFolder = f.getName().replace("_1.fq.gz", "")+"/";;//parentFolders[parentFolders.length-1]+"/";
				}
				writer.write("mkdir " + rootName + outputFolder+"\n");
				String fileName = "\""+f.getName()+"\"";
				fileName = fileName.replace("\\", "/");
				String fastqWithPath = f.getPath().replace("\\", "/");
				//String line = "kallisto quant --bias -t 4 -i " + indexFile + " -o " + rootName + outputFolder + " " + fileName + " " + fileName.replace(PairedStringForward, PairedStringReverse) + " &> "+ rootName + outputFolder+ "kallisto_"+outputFolder.replace("/","")+".err" +"\n";
				String line = "kallisto quant --bias -t 16 -i " + indexFile + " -o " + rootName + outputFolder + " \"" + fastqWithPath + "\" \"" + fastqWithPath.replace(PairedStringForward, PairedStringReverse) + "\" &> "+ rootName + outputFolder +outputFolder.replace("/","")+".err" +"\n";
				writer.write(line);
			}
			else if(f.getName().contains(".fq")) //single end
			{
				String outputFolder = f.getName().replace(".fq.gz","/");
				if(iris)
				{
					String[] parentFolders = f.getParent().split("\\\\");
					outputFolder = parentFolders[parentFolders.length-1]+"/";
				}
				writer.write("mkdir " +rootName + outputFolder+"\n");
				String fileName = "\""+f.getName()+"\"";
				fileName = fileName.replace("\\", "/");
				String fastqWithPath = f.getPath().replace("\\", "/");
//				if(fileName.contains("\'"))
//				{
//					fileName = fileName.replace("\'", "\\\'");
//					System.out.println(fileName);
//				}
				// --fragment-length=200 --sd=20 are standard values used. If more info is available about the fragment size use those
				//String line = "kallisto quant --bias -t 4 --single --fragment-length=200 --sd=20 -i " + indexFile + " -o " + rootName + outputFolder + " " + fileName + " &> "+ rootName + outputFolder+ "kallisto_"+outputFolder.replace("/","")+".err" +"\n";
				String line = "kallisto quant --bias -t 16 --single --fragment-length=200 --sd=20 -i " + indexFile + " -o " + rootName + outputFolder + " \"" + fastqWithPath + "\" &> "+ rootName + outputFolder+ "kallisto_"+outputFolder.replace("/","")+".err" +"\n";
				writer.write(line);
			}
		}
		writer.close();
		System.out.println("File written to: " + scriptName);
	}

	private static File[] getFileNames(String inputfilesFN) throws IOException 
	{
		
		BufferedReader reader= new BufferedReader(new FileReader(new File(inputfilesFN)));
		String line = null;
		int lines = 0;
		while(reader.readLine() != null)
		{
			lines++;
		}
		reader.close();
		reader= new BufferedReader(new FileReader(new File(inputfilesFN)));
		File[] files = new File[lines];
		int f = 0;
		while((line = reader.readLine()) != null)
		{
			files[f] = new File(line);
			f++;
		}
		reader.close();
		return files;
	}

}
