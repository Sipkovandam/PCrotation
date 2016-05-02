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
	public static void main(String[] args) throws IOException
	{
//		String folderFN = "E:/Groningen/Data/Monogenetic_disease_samples_RadBoud/";
//		String rootName = "/Volumes/Promise_RAID/GeneNetwork/Sipko/Monogenetic_disease_samples_RadBoud/";
//		String indexFile = "/Volumes/Promise_RAID/juha/PublicRNASeq/kallisto/data/hg19.v75.cdna.all.42.2.idx";
//		String scriptName = "E:/Groningen/Data/Monogenetic_disease_samples_RadBoud/KallistoRun.sh";
		
		String folderFN = "E:/Groningen/Data/Iris/allSamples_edited.txt";
		String rootName = "/local/groups/umcg-wijmenga/scr01/umcg-svandam/kallisto/Iris/";
		String indexFile = "/groups/umcg-wijmenga/scr01/umcg-svandam/kallisto/hg19.v75.cdna.all.42.2.idx";
		String scriptName = "E:/Groningen/Data/Iris/KallistoRunIrisSD200.sh";
		
		File[] files = null;
		if(folderFN.contains(".txt"))
			files = getFileNames(folderFN);
		else
			files = new File(folderFN).listFiles();
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(scriptName)));
		for(File f : files)
		{
			if(f.getName().contains("_R2_"))
				continue;
			else if(f.getName().contains("_R1_"))
			{
				String outputFolder = f.getName().split("_R1_")[0]+"/";
				writer.write("mkdir " + rootName + outputFolder+"\n");
				String fileName = "\""+f.getPath()+"\"";
				fileName = fileName.replace("\\", "/");
				String line = "kallisto quant --bias -t 4 -i " + indexFile + " -o " + rootName + outputFolder + " " + fileName + " " + fileName.replace("_R1_", "_R2_") + " &> "+ rootName + outputFolder + "mappingPercentages.txt" +"\n";
				writer.write(line);
			}
			else if(f.getName().contains(".fq")) //single end
			{
				String outputFolder = f.getName().replace(".fq.gz","/");
				writer.write("mkdir " +rootName + outputFolder+"\n");
				String fileName = "\""+f.getPath()+"\"";
				fileName = fileName.replace("\\", "/");
//				if(fileName.contains("\'"))
//				{
//					fileName = fileName.replace("\'", "\\\'");
//					System.out.println(fileName);
//				}
				// --fragment-length=200 --sd=20 are standard values used. If more info is available about the fragment size use those
				String line = "kallisto quant --bias -t 4 --single --fragment-length=200 --sd=200 -i " + indexFile + " -o " + rootName + outputFolder + " " + fileName + " &> "+ rootName + outputFolder + "mappingPercentages.txt" +"\n";
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
