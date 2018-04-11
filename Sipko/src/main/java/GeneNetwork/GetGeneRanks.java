package GeneNetwork;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;

import Tools.FileUtils;
import Tools.Script;

public class GetGeneRanks extends Script<GetGeneRanks>
{
	//folder containing the results from the Gene-network+GAVIN analysis
	String ranksCandidatesFolder = "";
	String ranksFolder = "E:/Groningen/Test/GeneNetwork.GetGeneRanks/BenchmarkSamples2/rankingAllgenes/";
	String sampleNameToCausalGeneFn = "";
	String writeFn = "";
	boolean blind = false;
	boolean spikeInCausalGene = true;
	
	@Override
	public void run()
	{
		try
		{
			File[] resultFns = new File(ranksCandidatesFolder).listFiles();
			
			HashMap<String,String> sampleName_To_CausalGene = FileUtils.readStringStringHash(sampleNameToCausalGeneFn);
			BufferedWriter resultsWriter = FileUtils.createWriter(writeFn);
			resultsWriter.write("Sample\tGene\tRank\tTotalGavinVariants\n");
			for(File resultFn : resultFns)
			{
				String sampleName = FileUtils.removeExtention(resultFn.getName());
				log("sampleName= " + sampleName);
				if(!sampleName_To_CausalGene.containsKey(sampleName))
					continue;
				
				String[] causalGenes = sampleName_To_CausalGene.get(sampleName).split("/");
				
				for(String causalGene:causalGenes)
				{
					log("SampleName= " + sampleName + "\tcausalGene = " + causalGene);
					
					String resultFnAllRanks= ranksFolder+resultFn.getName();
					int geneNetworkCausalGeneRank = FileUtils.getRankInFile(resultFnAllRanks, causalGene);
					
					causalGene=FileUtils.removeWhiteSpace(causalGene);
					String rank[] = getRankInFile(resultFn , causalGene, geneNetworkCausalGeneRank);
					
					//write the ranks
					if(blind)
					{
						log("sample= " + sampleName + "\trank= " + rank[0]+"/"+rank[1]+ "\t"+ rank[2]);
						resultsWriter.write(sampleName+"\t" + rank[0] + "\t" + rank[1] + "\t"+ rank[2] + "\n");
					}	
					else
					{
						log("sample= " + sampleName + "\tgene= " +causalGene + "\trank= " + rank[0]+"/"+rank[1]+ "\t"+ rank[2]);
						resultsWriter.write(sampleName+"\t"+causalGene+"\t" + rank[0] + "\t" + rank[1] + "\t"+ rank[2] + "\n");
					}
				}
			}
			resultsWriter.close();
		}catch(Exception e){e.printStackTrace();}
		log("done, file written to:\t" + writeFn);
	}

	private String[] getRankInFile(	File resultFn,
								String causalGene, int geneNetworkCausalGeneRank) throws FileNotFoundException, IOException
	{
		BufferedReader reader= FileUtils.createReader(resultFn.getAbsolutePath());
		String line = reader.readLine();
		
		String[] colHeaders=line.split("\t");
		
		int[] rank = new int[2];
		rank[0] =-1;//rank of causal gene
		rank[1] = 1;//total number of mutations
		String zScoresPerTerm="";
		while((line = reader.readLine())!=null)
		{
			if(line.contains(causalGene))
			{
				rank[0]=rank[1];
				
				//get HPO scores
				String[] columns = line.split("\t");

				zScoresPerTerm=addZscoresPerTerm(zScoresPerTerm, colHeaders, columns);
			}
			else if(rank[0]==-1 && spikeInCausalGene)
			{
				//get HPO scores
				String[] columns = line.split("\t");
				int geneNetworkRank= Integer.parseInt(columns[2]);
				if(geneNetworkCausalGeneRank<geneNetworkRank)
				{
					log("Gnrank= " + geneNetworkRank + "\tcausalGNrank=" + geneNetworkCausalGeneRank);
					rank[0]=rank[1];
					zScoresPerTerm=addZscoresPerTerm(zScoresPerTerm, colHeaders, columns);
				}
			}
			
			
			rank[1]++;
		}	
		reader.close();
		String[] results = new String[3];
		results[0]=Double.toString(rank[0]);
		results[1]=Double.toString(rank[1]);
		results[2]=zScoresPerTerm;
		return results;
	}

	private String addZscoresPerTerm(String zScoresPerTerm, String[] colHeaders, String[] columns)
	{
		zScoresPerTerm+=colHeaders[3]+"="+columns[3];
		for(int c=4; c< columns.length;c++)
		{
			zScoresPerTerm+="|"+colHeaders[c]+"="+columns[c];
		}
		return zScoresPerTerm;
	}
}
