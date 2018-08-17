package GeneNetwork;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.HashMap;

import Tools.FileUtils;
import Tools.Script;

public class UnannotatedSampleZscoreCounter extends Script<UnannotatedSampleZscoreCounter>
{
	//This script counts the number of times a gene has a z-score > X in a file.
	
	String inputFoldersFns="E:/Groningen/Papers/Gene Network Publication/Supplements/UnsolvedPatients/rankingCandidateGenesAnne/,E:/Groningen/Papers/Gene Network Publication/Supplements/UnsolvedPatients/PrioritisationsCardioEdgar/";
	String writeFn="";
//	String identifierToHpo="";
	double cutoff = 5;
	
	@Override
	public void run()
	{
		try
		{
//			HashMap<String,String> identifier_To_Hpo = FileUtils.readStringStringHash(identifierToHpo);
			String[] inputFolders = inputFoldersFns.split(",");
			BufferedWriter writer = FileUtils.createWriter(writeFn);
			String header = "Sample\tHPO terms\tn_Z_Above"+cutoff;
			writer.write(header+"\n");
			
			for(String inputFolder: inputFolders)
			{
				File[] fns= new File(inputFolder).listFiles();
				
				for(File fn : fns)
				{
					log("Fn=\t" + fn);
					BufferedReader zscoreReader = FileUtils.createReader(fn.getAbsolutePath());
					String sampleName= FileUtils.removeExtention(fn.getName());
					
					String line = zscoreReader.readLine();//header
					String[] eles = line.split("\t"); 
					String hpoTerms = getHpoTerms(eles);
					
					int n=0;
					while((line=zscoreReader.readLine())!=null)
					{
						double zScore=Double.parseDouble(line.split("\t")[3]);
						if(zScore < cutoff)
							break;
						n++;
					}
					String writeLine = sampleName +"\t"+ hpoTerms + "\t" + n;
					writer.write(writeLine+"\n");
					zscoreReader.close();
				}
			}
			writer.close();
			
		}catch(Exception e){e.printStackTrace();}
	}


	private String getHpoTerms(String[] eles)
	{
		String hpoTerms=eles[4];
		int length = eles.length;
		for(int i = 5; i<eles.length-1;i++)
		{
			hpoTerms+=","+eles[i];
		}
		return hpoTerms;
	}
}
