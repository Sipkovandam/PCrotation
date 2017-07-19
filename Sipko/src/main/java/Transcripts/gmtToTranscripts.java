package Transcripts;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import MatrixScripts.MatrixStruct;

public class gmtToTranscripts 
{
	//this script converts all the ensembl GENE IDs to ensembl Transcript IDs
	//1. if a gene has multiple transcripts only the coding transcripts are included into the resulting gmt file
	//If a gene only has non-coding transcripts (but more than 1) it is excluded. This may not necessarily be desirable... 
	
	static String ensgToEnstFN = "E:/Sipko/Thesis resubmission/Reactome Files/hg19.v75.cdna.all.enst2ensg.txt";
	static String gmtFN = "E:/Sipko/Thesis resubmission/ROC curve/ReactomePathways.ensg.gmt";
	static String writeFN = null;
	static String transcriptsToKeep = "E:/Sipko/Thesis resubmission/Reactome Files/ENST00000372179.txt";//optional
	static String cdsFN = "E:/Sipko/Thesis resubmission/Reactome Files/CodingSequences.txt";// (coding sequences obtained from biomart (col1:transcripts col2: cds)
	
	public static void main(String[] args) throws IOException
	{
		if(writeFN== null)
			writeFN = gmtFN.replace(".gmt", "_transcripts.gmt");
		
		if(transcriptsToKeep != null)
			ensgToEnstFN = keepGenes(transcriptsToKeep, ensgToEnstFN);
		
		createTranscriptGMT();
		
		
	}

	private static void createTranscriptGMT() throws IOException {
		Hashtable<String,ArrayList<String>> ensgToEnst = readConversionHash(ensgToEnstFN);
		BufferedReader reader = new BufferedReader(new FileReader(new File(gmtFN)));
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(writeFN)));
		Hashtable<String,ArrayList<String>> enstCoding = readConversionHash(cdsFN,0,1);
		System.out.println("enstCoding.length = " + enstCoding.size());
		
		String line = null;
		while((line=reader.readLine())!=null)
		{
			String[] eles = line.split("\t");//genes start from index 2
			String tsLine = eles[0]+"\t"+eles[1];
			int n = 0;
			for(int e = 2; e < eles.length; e++)
			{
				String gene = eles[e];
				ArrayList<String> transcripts = ensgToEnst.get(gene);
				if(transcripts == null)
					continue;
				for(int t = 0; t < transcripts.size(); t++)
				{
					String transcript = transcripts.get(t);
					if(transcripts.size()>1 && !enstCoding.containsKey(transcript))//if there are mutliple transcripts for this gene and there is no coding sequence for this transcript, continue
						continue;
					tsLine+="\t"+transcript;
					n++;
				}
			}
			if(n != 0)
				writer.write(tsLine+"\n");
		}
		reader.close();
		writer.close();
		System.out.println("Done! File written to:" + writeFN);
	}

	private static String keepGenes(String transcriptsToKeep2, String ensgToEnstFN2) throws IOException 
	{
		MatrixStruct genesToKeep = new MatrixStruct(transcriptsToKeep2);
		BufferedReader reader = new BufferedReader(new FileReader (new File(ensgToEnstFN2)));
		String writeFN = ensgToEnstFN2.replace(".txt", "_IDsInMap.txt");
		BufferedWriter writer = new BufferedWriter(new FileWriter (new File(writeFN)));
		
		String line = null;
		while((line = reader.readLine()) != null)
		{
			String[] eles = line.split("\t");
			//System.out.println(eles[0]);
			if(genesToKeep.rowHash.containsKey(eles[0]))
				writer.write(line+"\n");
		}
		
		writer.close();
		reader.close();
		System.out.println("File written to:" + writeFN);
		return writeFN;
	}

	private static Hashtable<String,ArrayList<String>> readConversionHash(String convertFile) throws IOException
	{
		return readConversionHash(convertFile,1,0);
	}
	
	private static Hashtable<String,ArrayList<String>> readConversionHash(String convertFile, int col1, int col2) throws IOException 
	{
		Hashtable<String,ArrayList<String>> ensgToEnst = new Hashtable<String,ArrayList<String>>();
		BufferedReader reader = new BufferedReader(new FileReader(new File(convertFile)));
		String line = null;
		while((line=reader.readLine())!= null)
		{
			String[] eles = line.split("\t");
			if(eles.length<2)
				continue;
			String gene = eles[col1];
			String transcript = eles[col2];
			ArrayList<String> transcripts = ensgToEnst.get(gene);
			if(transcripts == null)
				transcripts = new ArrayList<String>();
			
			transcripts.add(transcript);
			ensgToEnst.put(gene, transcripts);
		}
		reader.close();
		return ensgToEnst;
	}
}
