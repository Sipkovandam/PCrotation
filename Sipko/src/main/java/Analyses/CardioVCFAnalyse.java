package Analyses;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;

public class CardioVCFAnalyse 
{
	//gets all variants that exist in at least 1 individual (individuals start from column 9)
	//only variants with HIGH or MODERATE impact included
	//only variants that PASS
	//AFcutoff is the Allele frequency cutoff in the 1000.genomes project)
	static double AFcutoff = 0.01;
	
	public static void main(String[] args) throws NumberFormatException, IOException
	{
		String fileName = "E:/Groningen/Data/Tessa/Golden_Standard_WES_Cardio.final.vcf";
		String writeFN = "E:/Groningen/Data/Tessa/Golden_Standard_WES_Cardio.final.RareVariants"+AFcutoff+"_High.txt";
		File file = new File(fileName);
				
		BufferedReader reader = new BufferedReader(new FileReader(file));
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(writeFN)));
		String line = null;
//		String line = reader.readLine();
//		String[] headers = line.split("\t");
//		Hashtable<String,Integer> headerToCol = new Hashtable<String,Integer>();
//		for(int i=0; i < headers.length; i++)
//			headerToCol.put(headers[i], i);
		String header = "ensembl_geneID\tTranscriptName\tgeneName\timpact\talleleFrequency\tHomozygoot\tHeterozygoot\tBiotype";
		System.out.println(header);
		writer.write(header);
		
		while((line=reader.readLine())!=null)
		{
			String[] cells = line.split("\t");
			if(cells.length<7 || !cells[6].contains("PASS"))
				continue;
			
//			int col = headerToCol.get("");
			//System.out.println(cells[7]);
			String impact = get("SNPEFF_IMPACT=", cells);
			
			if(impact ==null || (!impact.contains("HIGH")))
				continue;
			String alleleFrequency = get("dbNSFP_1000Gp1_EUR_AF=", cells);
			if(alleleFrequency == null || alleleFrequency.equals(".") || alleleFrequency.equals(",.,"))
				continue;
			
			
			try
			{				
				if(Double.parseDouble(alleleFrequency) > AFcutoff)
					continue;
			}catch(Exception e)
			{
//				e.printStackTrace();
				continue;
			}
			
			String transcriptName = get("SNPEFF_TRANSCRIPT_ID=", cells);
			String geneName = get("SNPEFF_GENE_NAME=", cells);
			String biotype = get("SNPEFF_GENE_BIOTYPE=", cells);
			String ensembl_geneID = get("dbNSFP_Ensembl_geneid=", cells);
			
//			if(transcriptName == null || geneName == null || biotype == null)
//			{
//				System.out.println(transcriptName+"\t"+geneName+"\t"+impact+"\t"+alleleFrequency+"\t"+biotype);
//				continue;
//			}
			
			int patients = cells.length;
			int startCol = 9;
			int homozygoot = 0;
			int heterozygoot = 0;
			for(int p=0+startCol; p< patients; p++)
			{
				if(cells[p].contains("1/1"))
					homozygoot++;
				else if(cells[p].contains("/1") || cells[p].contains("1/"))
					heterozygoot++;
			}
			String outputLine = ensembl_geneID+"\t"+transcriptName+"\t"+geneName+"\t"+impact+"\t"+alleleFrequency+"\t"+homozygoot+"\t"+heterozygoot+"\t"+biotype;
			System.out.println(outputLine);
			writer.write(outputLine+"\n");
		}	
		reader.close();
		writer.close();
	}

	private static String get(String seachString, String[] cells) {
		if(!cells[7].contains(seachString))
			return null;
		return cells[7].split(seachString)[1].split(";")[0];
	}
}
