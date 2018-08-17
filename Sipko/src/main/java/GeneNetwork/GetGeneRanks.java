package GeneNetwork;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import MatrixScripts.GetRows;
import Tools.FileUtils;
import Tools.Script;

public class GetGeneRanks extends Script<GetGeneRanks>
{
	//folder containing the results from the Gene-network+GAVIN analysis
	String ranksCandidatesFolder = "";
	String ranksFolder = "E:/Groningen/Test/GeneNetwork.GetGeneRanks/BenchmarkSamples2/rankingAllgenes/";
	String sampleNameToCausalGeneFn = "";
	String writeFn = "";
	String diseaseGenesFn= "E:/Groningen/Data/Annotation/GRCh37/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt_genesInPathways.txt";
	boolean blind = false;
	boolean spikeInCausalGene = true;
	String dnaNumberToRvcf = "E:/Groningen/Test/GeneNetwork.GetGeneRanks/BenchmarkSamples2/dnaNumbersToRvcf_SipkoComputer.txt";
	
	@Override
	public void run()
	{
		try
		{
			if(diseaseGenesFn!=null)
				getDiseaseGenesFromFiles(ranksCandidatesFolder);
			
//			if(dnaNumberToRvcf!=null)
//				addVariantsToFiles(ranksCandidatesFolder, dnaNumberToRvcf);
			
			getRanksFromFolder(ranksCandidatesFolder, writeFn, null);
			
		}catch(Exception e){e.printStackTrace();}
		log("done, file written to:\t" + writeFn);
	}

	private void addVariantsToFiles(String ranksCandidatesFolder, String dnaNumberToRvcf) throws FileNotFoundException, IOException
	{
		String addToFileName = "_VariantsAdded";
		String writeFolderName = FileUtils.makeFolderNameEndWithSlash(ranksCandidatesFolder)+"VariantsAdded/";
		FileUtils.makeDir(writeFolderName);
		File folder = new File(ranksCandidatesFolder);
		
		HashMap<String, String> dnaNumber_To_GavinFn = FileUtils.readStringStringHash(dnaNumberToRvcf);
		
		for(File file : folder.listFiles())
		{
			if(file.isDirectory())
				continue;
			String fn = file.getName();
			String dnaNumber = FileUtils.removeExtention(fn);
			
			String writeFn = writeFolderName+FileUtils.removeExtention(fn)+addToFileName+".txt";
			String gavinFn = dnaNumber_To_GavinFn.get(dnaNumber);
			log("DnaNumber = \t" + dnaNumber + "\t" +gavinFn);
			HashMap<String,ArrayList<String>> gene_To_GavinLines = parseRvcf(gavinFn);
			
			BufferedReader reader = FileUtils.createReader(file.getAbsolutePath());
			BufferedWriter writer = FileUtils.createWriter(writeFn);
			//write header
			String line = reader.readLine();
			writer.write(line+"\n");
			
			while((line=reader.readLine())!=null)
			{
				String[] eles = line.split("\t");
				String geneSymbol = eles[1];
				ArrayList<String> gavinLines = gene_To_GavinLines.get(geneSymbol);
				if(gavinLines==null)//
				{
					log("Gavin line is missing for gene:\t"+geneSymbol+"\t this should not be possible! Exiting");
					continue;
					//System.exit(1);
				}
				
				StringBuilder newLineBuilder = new StringBuilder();
				newLineBuilder.append(line);
				for(String gavinLine: gavinLines)
				{
					newLineBuilder.append("\t");
					newLineBuilder.append(gavinLine);
				}
				
				newLineBuilder.append("\n");
				writer.write(newLineBuilder.toString());
			}
			

			writer.close();
			reader.close();
		}
		
	}

	private HashMap<String, ArrayList<String>> parseRvcf(String gavinFn) throws FileNotFoundException, IOException
	{
		HashMap<String, ArrayList<String>> gene_To_GavinLines = new HashMap<String, ArrayList<String>>();
		BufferedReader gavinReader = FileUtils.createReader(gavinFn);
		String header = gavinReader.readLine();
		String line = null;
		while((line=gavinReader.readLine())!=null)
		{
			String geneSymbol = line.split("\t")[8];
			ArrayList<String> geneLines = gene_To_GavinLines.get(geneSymbol);
			if(geneLines==null)
				geneLines=new ArrayList<String>();
			geneLines.add(line);
			gene_To_GavinLines.put(geneSymbol, geneLines);
		}
		
		return gene_To_GavinLines;
	}

	private void getRanksFromFolder(String ranksCandidatesFolder, String writeFn, String addedToFileName) throws FileNotFoundException, IOException
	{
		File[] resultFns = new File(ranksCandidatesFolder).listFiles();
		
		HashMap<String,String> sampleName_To_CausalGene = FileUtils.readStringStringHash(sampleNameToCausalGeneFn);
		BufferedWriter resultsWriter = FileUtils.createWriter(writeFn);
		resultsWriter.write("Sample\tGene\tRank\tTotalGavinVariants\n");
		for(File resultFn : resultFns)
		{
			String sampleName = FileUtils.removeExtention(resultFn.getName());

			if(addedToFileName!=null)//beetje brak, maar wel snel :S
				sampleName=sampleName.replace(addedToFileName, "");
			
			log("sampleName= " + sampleName);
			if(!sampleName_To_CausalGene.containsKey(sampleName))
				continue;
			
			String[] causalGenes = sampleName_To_CausalGene.get(sampleName).split("/");
			
			for(String causalGene:causalGenes)
			{
				log("SampleName= " + sampleName + "\tcausalGene = " + causalGene);
				
				String resultFnAllRanks= ranksFolder+resultFn.getName();
				if(addedToFileName!=null)
					resultFnAllRanks=resultFnAllRanks.replace(addedToFileName, "");
				
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
	}

	private void getDiseaseGenesFromFiles(String ranksCandidatesFolder) throws FileNotFoundException, IOException
	{
		String addToFileName = "_diseaseGenes";
		String writeFolderName = FileUtils.makeFolderNameEndWithSlash(ranksCandidatesFolder)+"diseaseGenesOnly/";
		FileUtils.makeDir(writeFolderName);
		File folder = new File(ranksCandidatesFolder);
		GetRows rowGetter = new GetRows();
		rowGetter.setRowsToGetFn(diseaseGenesFn);
		
		for(File file : folder.listFiles())
		{
			if(file.isDirectory())
				continue;
			String writeFn = writeFolderName+FileUtils.removeExtention(file.getName())+addToFileName+".txt";
			rowGetter.setFileName(file.getAbsolutePath());
			rowGetter.setWriteName(writeFn);
			rowGetter.setOriginalOrder(true);
			rowGetter.run();
		}
		
		getRanksFromFolder(writeFolderName, FileUtils.removeExtention(writeFn)+"_diseaseGenesOnly.txt", addToFileName);
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
			if(line.split("\t")[1].equals(causalGene))
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
		for(int c=4; c< 5;c++)
		{
			zScoresPerTerm+="|"+colHeaders[c]+"="+columns[c];
		}
		return zScoresPerTerm;
	}
}
