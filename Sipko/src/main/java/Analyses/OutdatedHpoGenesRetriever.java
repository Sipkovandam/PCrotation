package Analyses;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import MatrixScripts.MyMatrix;
import Tools.FileUtils;
import Tools.Script;

public class OutdatedHpoGenesRetriever extends Script<OutdatedHpoGenesRetriever>
{
	String omimFn = "E:/Groningen/Data/Annotation/GRCh37/OMIM/morbidmap2.txt";
	String hpoFn = "E:/Groningen/Data/Annotation/GRCh37/OMIM/ALL_SOURCES_FREQUENT_FEATURES_diseases_to_genes_to_phenotypes_OMIMonly.txt";
	String hpoZscoreFn = "E:/Groningen/Data/Juha/31-01-2018/HPO/hpo_predictions_bonSigOnly.txt.gz";
	String skewnessFn = "E:/Groningen/Students/Patrick/GeneNetwork/PathwayPredictions/hpo_predictions_genePredictability.txt";
	String ensembl_To_GeneSymbolFn = "E:/Groningen/Data/Annotation/GRCh37/ENSGToGeneNameV75.txt"; 
	String writeFnMissing = "E:/Groningen/Data/Annotation/GRCh37/OMIM/ALL_SOURCES_FREQUENT_FEATURES_diseases_to_genes_to_phenotypes_OMIMonly_MissingFromOmim.txt";
	String writeFnPresent= "E:/Groningen/Data/Annotation/GRCh37/OMIM/ALL_SOURCES_FREQUENT_FEATURES_diseases_to_genes_to_phenotypes_OMIMonly_PresentInOmim.txt";
	
	transient HashMap<String,HashSet<String>> gene_To_HpoTerms = null;
	transient HashMap<String, String> ensebml_To_GeneSymbol = null;
	transient HashMap<String, ArrayList<String>> geneSymbol_To_Ensebml = null;
	transient HashMap<String, String> geneSymbol_To_GeneSignificanceHpo = null;
	transient HashSet<String> omimIdsToIgnore = new HashSet<String>();//Basically: (4) the disorder is a chromosome deletion or duplication syndrome.
	transient HashMap<String,HashSet<String>> omim_To_PresentGenes = new HashMap<String,HashSet<String>>();
	transient HashMap<String,HashSet<String>> omim_To_MissingGenes = new HashMap<String,HashSet<String>>();
	transient MyMatrix hpoToZscores = null;
	
	@Override
	public void run()
	{
		try
		{
			init();
			
			//create a HashMaps of all genes with each omim term
			HashMap<String,HashMap<String,Integer>> omimTerm_To_genes = readOmimFile();
			
			//get all genes that are present in HPO that are NOT there in OMIM
			addMissingGenesToHash(omimTerm_To_genes, omim_To_PresentGenes, omim_To_MissingGenes);
			
			//write results
			writeMissingOmimToMissingGenes(omim_To_MissingGenes,writeFnMissing);
			
			//header = Omim	HP:000x	z-score	predictabilityScore
			writeResultTable(omim_To_MissingGenes,this.writeFnMissing);
			writeResultTable(omim_To_PresentGenes,this.writeFnPresent);
			
		}catch(Exception e){e.printStackTrace();}
		
		log("Finished, file written to: " + writeFnMissing);
	}

	private void writeResultTable(HashMap<String, HashSet<String>> omim_To_MissingGenes, String writeFn) throws FileNotFoundException, IOException
	{
		writeFn = FileUtils.removeExtention(writeFn)+"_allInfoTable.txt";
		BufferedWriter geneWriter = FileUtils.createWriter(writeFn);
		geneWriter.write("OmimNumber\tHpoTerm\tGeneSymbol\tEnsemblId\tZscore\tSignificanceScore\n");
		for(String omimTerm: omim_To_MissingGenes.keySet())
		{
			log("omimTerm=\t" + omimTerm);
			if(omimIdsToIgnore.contains(omimTerm))
				continue;
			
			HashSet<String> genes = omim_To_MissingGenes.get(omimTerm);
			for(String gene:genes)
			{				
				log("gene=\t" + gene);
				HashSet<String> hpoTerms = gene_To_HpoTerms.get(gene);				
				for(String hpoTerm : hpoTerms)
				{		
					double zScore = Double.NaN;
					//log("hpoterm=\t" + hpoTerm + "\t" + hpoToZscores.getRowHash().keys().nextElement() + "\t" + hpoToZscores.getColHash().keys().nextElement());
					if(hpoTerm.equals("HP:0006297"))
						log("hpoterm=\t" + hpoTerm + "\tGS=" + geneSymbol_To_Ensebml.containsKey(gene) + "\t" + hpoToZscores.getColHash().containsKey(hpoTerm));
					
//					if(geneSymbol_To_Ensebml.containsKey(gene))
//						log("geneEnsembl=\t" + geneSymbol_To_Ensebml.get(gene) + "\thasGene=\t" + hpoToZscores.getRowHash().containsKey(geneSymbol_To_Ensebml.get(gene)));
					
					ArrayList<String> ensemblIds= geneSymbol_To_Ensebml.get(gene);
					
					if(ensemblIds==null)
						continue;
					
					for(String ensemblId: ensemblIds)
					{
						if(hpoToZscores.getRowHash().containsKey(ensemblId) && hpoToZscores.getColHash().containsKey(hpoTerm))
						{	
							int row = hpoToZscores.getRowHash().get(ensemblId);
							int col = hpoToZscores.getColHash().get(hpoTerm);
							zScore = hpoToZscores.values[row][col];
							if(hpoTerm.equals("HP:0006297"))
								log("z-score=\t" + zScore);
						}
						StringBuilder lineBuilder = new StringBuilder();
						lineBuilder.append(omimTerm);
						lineBuilder.append("\t");
						lineBuilder.append(hpoTerm);
						lineBuilder.append("\t");
						lineBuilder.append(gene);
						lineBuilder.append("\t");
						lineBuilder.append(ensemblId);
						lineBuilder.append("\t");
						lineBuilder.append(zScore);
						lineBuilder.append("\t");
						
						//get gene significance
						String geneSignificance = "NaN";
						if(geneSymbol_To_GeneSignificanceHpo.containsKey(gene))
						{
							geneSignificance = geneSymbol_To_GeneSignificanceHpo.get(gene);
						}
						lineBuilder.append(geneSignificance);
						lineBuilder.append("\n");
						geneWriter.write(lineBuilder.toString());
					}
				}
			}
			
		}
		geneWriter.close();
		
		log("AllInfoTable table written to:\t" + writeFn);
	}

	private void writeMissingOmimToMissingGenes(HashMap<String, HashSet<String>> omim_To_MissingGenes,
												String writeFn) throws FileNotFoundException, IOException
	{
		BufferedWriter missingGeneWriter = FileUtils.createWriter(writeFn);
		missingGeneWriter.write("OmimNumber\tgenes\n");
		for(String omimTerm: omim_To_MissingGenes.keySet())
		{
			if(omimIdsToIgnore.contains(omimTerm))
				continue;
			StringBuilder lineBuilder = new StringBuilder();
			lineBuilder.append(omimTerm);
			
			HashSet<String> genes = omim_To_MissingGenes.get(omimTerm);
			for(String gene:genes)
			{
				lineBuilder.append("\t");
				lineBuilder.append(gene);
			}
			lineBuilder.append("\n");
			missingGeneWriter.write(lineBuilder.toString());
		}
		missingGeneWriter.close();
	}

	private void init() throws FileNotFoundException, IOException
	{
		if(writeFnMissing==null)
			writeFnMissing=FileUtils.removeExtention(hpoFn)+"_MissingFromOmim.txt";
		
		log("Reading in genes To HPO terms");
		gene_To_HpoTerms = readGene_To_HpoTerms(hpoFn);
		
		log("Reading in ensembl to gene symbol file");
		ensebml_To_GeneSymbol= FileUtils.readStringStringHash(ensembl_To_GeneSymbolFn);
		
		
		geneSymbol_To_Ensebml= FileUtils.readStringMultiStringArrayList(ensembl_To_GeneSymbolFn,1,0, false);
		
		log("Reading in geneSymbol to geneSignificanceHpo file");
		HashMap<String, String> ensebmlId_To_GeneSignificanceHpo = FileUtils.readStringStringHash(skewnessFn,0,2);
		geneSymbol_To_GeneSignificanceHpo= new HashMap<String, String>();
		log("Converting ensemblIDs to gene symbols");
		Set<String> genes = ensebmlId_To_GeneSignificanceHpo.keySet();
		for(String gene: genes)
		{
			if(ensebml_To_GeneSymbol.containsKey(gene))
			{
				String geneSymbol = ensebml_To_GeneSymbol.get(gene);
				geneSymbol_To_GeneSignificanceHpo.put(geneSymbol, ensebmlId_To_GeneSignificanceHpo.get(gene));
			}
		}
		
		log("Reading in z-score matrix");
		hpoToZscores = new MyMatrix(hpoZscoreFn);
		
	}

	private HashMap<String, HashSet<String>> readGene_To_HpoTerms(String hpoFn) throws FileNotFoundException, IOException
	{
		HashMap<String, HashSet<String>> gene_To_HpoTerm = new HashMap<String, HashSet<String>>();
		BufferedReader hpoReader = FileUtils.createReader(hpoFn);
		String line= hpoReader.readLine();//file has header
		while((line=hpoReader.readLine())!=null)
		{
			String[] cells = line.split("\t");
			String gene = cells[1];
			String hpoNumber = cells[3];
			gene_To_HpoTerm = FileUtils.addStringToHashSetInHash(gene,hpoNumber , gene_To_HpoTerm);
		}
		
		return gene_To_HpoTerm;
	}

	private HashMap<String, HashMap<String, Integer>> readOmimFile() throws FileNotFoundException, IOException
	{
		HashMap<String,HashMap<String,Integer>> omimTerm_To_genes = new HashMap<String,HashMap<String,Integer>>();
		BufferedReader omimReader=FileUtils.createReader(omimFn);
		
		String line = omimReader.readLine();//header line
		
		while((line=omimReader.readLine())!=null)
		{
			String[] cells = line.replace("\"", "").split("\t");
			String[] omimSplit = cells[0].split(" ");
			String omimNumber = omimSplit[omimSplit.length-2];
			String omimCode = omimSplit[omimSplit.length-1];//(4) means the disorder is a chromosome deletion or duplication syndrome
			//log("omimnumber=\t" + omimNumber);
			String[] genes= cells[1].split(",");
			
			HashMap<String,Integer> geneHash = omimTerm_To_genes.get(omimNumber);
			
			//log("OmimCode =" + omimCode);
			if(omimCode.endsWith("(4)"))
			{
				omimIdsToIgnore.add(omimNumber);
			}
			
			if(geneHash==null)
				geneHash=new HashMap<String,Integer>();
			
			for(String gene: genes)
			{
				if(gene.equals("NPR2L"))
					log("NPR2L:\t"+ omimNumber);
				geneHash.put(gene, 1);
			}
			
			omimTerm_To_genes.put(omimNumber, geneHash);
			
			
		}
		omimReader.close();
		return omimTerm_To_genes;
	}

	private void addMissingGenesToHash(HashMap<String, HashMap<String, Integer>> omimTerm_To_genes, HashMap<String, HashSet<String>> omim_To_PresentGenes, HashMap<String, HashSet<String>> omim_To_MissingGenes) throws FileNotFoundException, IOException
	{
		BufferedReader hpoReader = FileUtils.createReader(hpoFn);
		String line = hpoReader.readLine();//header
		while((line=hpoReader.readLine())!=null)
		{
			String[] cells = line.split("\t");
			String omimNumber=cells[0].replace("OMIM:", "");
			
			String gene = cells[1];
			
			if(!omimTerm_To_genes.containsKey(omimNumber))
				continue;
			
			boolean isPresent = omimTerm_To_genes.get(omimNumber).containsKey(gene);

			//just skip these genes and consider them as OK
			if(!omimTerm_To_genes.containsKey(omimNumber))
				isPresent=true;
			
			log("omimnumber=\t" + omimNumber + "\tgene=\t" + gene + "\t" + isPresent);
			
			if(!isPresent)
				omim_To_MissingGenes= FileUtils.addStringToHashSetInHash(omimNumber, gene, omim_To_MissingGenes);
			else
				omim_To_PresentGenes= FileUtils.addStringToHashSetInHash(omimNumber, gene, omim_To_PresentGenes);
			
		}
		hpoReader.close();
	}
}
