package GeneNetwork;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;

import MatrixScripts.MyMatrix;
import Tools.FileUtils;
import Tools.Script;

public class ExomizerRankRetriever extends Script<ExomizerRankRetriever>
{
	String exomizerOutputFolder = "/groups/umcg-gdio/tmp03/umcg-svandam/Data/VCF/Mixed/ExomizerSettings/"; 
	String writeFn = "/groups/umcg-gdio/tmp03/umcg-svandam/Data/VCF/Mixed/ExomizerRanks.txt";
	String dna_To_CausalGeneFn = "/groups/umcg-gdio/tmp03/umcg-svandam/Data/VCF/Mixed/Dna_To_CausalGene.txt";
	
	String gnRanksFolder = "E:/Groningen/Test/GeneNetwork.GetGeneRanks/BenchmarkSamples2/rankingAllgenes/";
	
	@Override
	public void run()
	{
		try
		{
			File exomizerFolder = new File(exomizerOutputFolder);
			HashMap<String,String> dna_To_CausalGene = FileUtils.readStringStringHash(dna_To_CausalGeneFn);
			BufferedWriter outputWriter = FileUtils.createWriter(writeFn);
			outputWriter.write("DNA number\tGene\tRank\tTotalListed\tgeneNetworkRank\n");
			
			BufferedWriter failedWriter = FileUtils.createWriter(writeFn.replace(".txt", "_failed.txt"));
			for(File folder : exomizerFolder.listFiles())
			{
				if(!folder.isDirectory())
					continue;
				
				String fn = folder.getAbsolutePath()+"/combined_2Columns.txt";
				
				log("fn=" + fn);
				
				MyMatrix exomizerResult = new MyMatrix(fn);
				
				int rank = -1;
				
				if(exomizerResult.values[0][0]==exomizerResult.values[exomizerResult.rows()-1][0] || new File(fn).length()==0)
				{
					failedWriter.write(fn+"\n");
				}
				
				exomizerResult.sortCol(0);
				
				String dnaNumber=folder.getName();
				
				log("DNA number =\t" + dnaNumber);
				String[] causalGenes = dna_To_CausalGene.get(dnaNumber).split("/");
				
				for(String causalGene:causalGenes)
				{
					int geneNetworkRank = getGeneNetworkRankAmongExomiserGenes(causalGene, exomizerResult.getRowHash(false), dnaNumber);
					
					log("CausalGene =\t" + causalGene);
				
					if(causalGene==null)
						log("Warning! DnaNumber not found in causalGeneFile:\t" + dnaNumber);		
					
					if(rank==-1 && exomizerResult.getRowHash(false).containsKey(causalGene))
					{
						rank = correctForEqualRanks(exomizerResult, causalGene);
					}
					
					String outputLine = dnaNumber+"\t"+causalGene+"\t"+ rank +"\t" + exomizerResult.rows()+"\t"+ geneNetworkRank +"\n";
					
					log("Rank=\t" + outputLine);
					outputWriter.write(outputLine);
					
					
				}
			}
			failedWriter.close();
			outputWriter.close();
			log("Finished, file written at:\t" + writeFn);
		}catch(Exception e){e.printStackTrace();}
	}

	private int correctForEqualRanks(MyMatrix exomizerResult, String causalGene)
	{
		int rank = exomizerResult.getRowHash().get(causalGene);

		double value=exomizerResult.values[rank][0];
		
		int i = rank;
		while(i >= 0 && exomizerResult.values[i][0]==value)
		{
			i--;
		}
		int j = rank;
		while(exomizerResult.values[j][0]==value && j < exomizerResult.values[j].length)
		{
			
			j++;
		}
		int avgRank=(i+j)/2+1;
		
		return avgRank;
	}

	private int getGeneNetworkRankAmongExomiserGenes(String causalGene, Hashtable<String, Integer> rowHash, String dnaNumber) throws FileNotFoundException, IOException
	{
		int gnRank=1;
		String gnFn=FileUtils.makeFolderNameEndWithSlash(gnRanksFolder)+ dnaNumber + ".txt";
		log("Retrieving gn rank from:" + gnFn);
		BufferedReader gnReader = FileUtils.createReader(gnFn);
		String line = gnReader.readLine();
		while((line = gnReader.readLine())!=null)
		{
			String gene = line.split("\t")[1];

			if(gene.equals(causalGene))
			{
				log("found!");
				gnReader.close();
				return gnRank;
			}
			if(rowHash.containsKey(gene))
				gnRank++;
		}
		
		gnReader.close();
		return -1;
	}

	
	
}
