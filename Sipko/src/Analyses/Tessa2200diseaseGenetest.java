package Analyses;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Hashtable;

import pca.MatrixStruct;

public class Tessa2200diseaseGenetest 
{
	//This script can be used to test the quality of the network by determining its ability to identify known disease genes.
	//This script determines where known disease genes rank based on their OMIM terms (which correspond to a number of HPO terms)
	
	static String zScoresGenesVSHPOFN = null;//zScoreFN
	static String diseaseGeneFN = null;//has the disease genes (gene symbols) and corresponding HPO terms (dp)
	static String geneSymbolToEnsemblFN = null; //idToEnsg
	static String clinvar = null;//clinvar
	
//	static String diseaseGeneFN = "E:/Groningen/Scripts/Tessa/files/HPO_ALL_SOURCES_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes.txt";//has the disease genes (gene symbols) and corresponding HPO terms
//	static String geneSymbolToEnsemblFN = "E:/Groningen/Data/PublicSamples/ENSGToGeneNameHGNCBiotypeChromosomeStartStopStrandAndDescriptionV75.txt.filtered.txt";
//	static String clinvar = "E:/Groningen/Scripts/Tessa/files/ClinVar_variant_summary.txt";
	
//	static String zScoresGenesVSHPOFN = "E:/Groningen/Data/PublicSamples/06-2016/5Kmatrix/Predictions/TermGeneZScores-2016-06-11HPO-unscaled.txt";
//	static String zScoresGenesVSHPOFN = "E:/Groningen/Data/PublicSamples/06-2016/22214matrix_reproduce/Prediction/Predictions/TermGeneZScores-2016-06-18HPO-unscaled.txt";
//	static String zScoresGenesVSHPOFN = "E:/Groningen/Data/PublicSamples/06-2016/5Kmatrix_reproduce/Prediction/Predictions/TermGeneZScores-2016-06-17HPO-unscaled.txt";
//	static String zScoresGenesVSHPOFN = "E:/Groningen/Data/PublicSamples/06-2016/Juha_QuantNorm_Covariance/Prediction/Predictions/TermGeneZScores-2016-06-10HPO-unscaled.txt";
//	static String zScoresGenesVSHPOFN = "E:/Groningen/Data/PublicSamples/06-2016/Juha_QuantNorm_Cov_zTransformed/Predictions/TermGeneZScores-2016-06-15HPO-unscaled.txt";
//	static String zScoresGenesVSHPOFN = "E:/Groningen/Data/PublicSamples/05-2016/Tessa/est_counts_nocancernocellline_Rlog_covariance_Chr1-22_add0.0BeforeGeoMean/Prediction/Predictions/TermGeneZScores-2016-06-09HPO-unscaled.txt";
//	static String zScoresGenesVSHPOFN = "E:/Groningen/Data/PublicSamples/06-2016/22214_Quantnorm_Cov_GCcorrected/Prediction/Predictions/TermGeneZScores-2016-06-20HPO-unscaled.txt";
	
//	E:\Groningen\Data\PublicSamples\05-2016\Tessa\22214Samples_Rlog_correlation_+1log\Prediction\Predicitions\
//	E:/Groningen/Data/PublicSamples/05-2016/Tessa/22214Samples_Rlog_correlation_+1log/Prediction/Predicitions/TermGeneZScores-2016-06-22HPO-unscaled.txt
//	E:/Groningen/Data/PublicSamples/06-2016/22214Samples_Rlog_covariance/Prediction/Predictions/TermGeneZScores-2016-06-22HPO-unscaled.txt
	
	//static String zScoresGenesVSHPOFN = "E:/Groningen/Data/PublicSamples/06-2016/22214matrix_reproduce/Prediction/Predictions/TermGeneZScores-2016-06-18HPO-unscaled_HPOshopping.txt";
	
	static int maxHPOsPerGene = -1;// -1 means this argument is ignored
	static int maxGenesPerOmim = -1;// -1 means this argument is ignored
	
	public static void main(String[] args) throws Exception
	{
		checkArgs(args);
		
		String writeFN = zScoresGenesVSHPOFN.replace(".txt", "_PredictionRankings.txt");
		String IDS_ConvertedFN = convertGStoEng(1,diseaseGeneFN);//first argument is the column
		String clinvarFN = convertGStoEng(4,clinvar);//first argument is the column
		
		//make a hashtable that contains all the Disease genes and the corresponding HPOs
		Hashtable<String, ArrayList<String>> OMIM_HPOs = readOMIM_HPOs(IDS_ConvertedFN);
		System.out.println("Number of OMIM disease terms = " + OMIM_HPOs.size());
		
		//make a hash with gene and corresponding OMIM
		Hashtable<String, HashSet<String>> geneToOMIM = new Hashtable<String, HashSet<String>>();
		Hashtable<String, Integer> OMIMgenes = getOMIMgenes(clinvarFN,geneToOMIM);
		
		//for each disease gene determine the rank per OMIM term and write the result
		MatrixStruct zScores = new MatrixStruct(zScoresGenesVSHPOFN);//HPO terms are on the rows, genes on the columns
		Enumeration<String> geneList = geneToOMIM.keys();
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(writeFN)));
		writer.write("GeneName\tRank\tOMIM\tOMIMgenes\tHPOs\n");
		int n = 0;
		int rankTop1000 = 0;
		while(geneList.hasMoreElements())
		{
			String gene = geneList.nextElement();
			HashSet<String> OMIMS = geneToOMIM.get(gene);
			for(String OMIM : OMIMS)
			{
				//System.out.println("number of genes to omimterm = " + OMIMgenes.get(OMIM));//OMIMgenes
				if(maxGenesPerOmim > 0 && OMIMgenes.get(OMIM) > maxGenesPerOmim)
					continue;
				ArrayList<String> HPOs = OMIM_HPOs.get(OMIM);
				if(HPOs == null || (maxHPOsPerGene > 0 && HPOs.size() > maxHPOsPerGene))
					continue;
				int rank = getRank(gene,HPOs,zScores);
				String writeLine = gene+"\t"+ rank +"\t"+ OMIM + "\t" + OMIMgenes.get(OMIM) + "\t" + HPOs.size() + "\n";
				if(rank != -1)
				{
					if(rank < 1001)
						rankTop1000++;
					System.out.println("n= " + n +"\tGene\t" + gene + "\tOMIM\t"+ OMIM + "\tRank\t" + rank + "\trankTop1000=" + rankTop1000);
				}
				writer.write(writeLine);
				n++;
			}
		}
		
		writer.close();
		System.out.println("File written to " + writeFN);
	}

	private static Hashtable<String, Integer> getOMIMgenes(String clinvarFN, Hashtable<String, HashSet<String>> geneToOMIM) throws IOException {
		Hashtable<String, Integer> OMIMgenes = new Hashtable<String, Integer>();
		BufferedReader geneReader = new BufferedReader(new FileReader(new File(clinvarFN)));
		String line = geneReader.readLine();
		
		//for each gene get the list of omims that belong to it
		while((line=geneReader.readLine()) != null)
		{
			String eles[] = line.split("\t");
			String geneName = eles[4];
			String[] OMIMeles = eles[10].split(",|;");
			for(String ele : OMIMeles)
			{
				if(!ele.contains("OMIM"))//if there is no omim ID on this line
					continue;
				//System.out.println(geneName + "\t" +ele);
				HashSet<String> geneOMIMS = geneToOMIM.get(geneName);
				if(geneOMIMS == null)
					geneOMIMS = new HashSet<String>();
				geneOMIMS.add(ele);
				geneToOMIM.put(geneName, geneOMIMS);
				if(!OMIMgenes.containsKey(ele))
					OMIMgenes.put(ele, 1);
				else
					OMIMgenes.put(ele, OMIMgenes.get(ele)+1);
			}
		}
		System.out.println("Number of disease genes = " + geneToOMIM.size());
		geneReader.close();
		return OMIMgenes;
	}

	private static int getRank(String gene, ArrayList<String> HPOs, MatrixStruct zScores) 
	{
		int rank = -1;
		if(!zScores.colHash.containsKey(gene))
			return -1;
		MatrixStruct geneZscores = null;
		for(int h = 0; h < HPOs.size(); h++)
		{
			if(!zScores.rowHash.containsKey(HPOs.get(h)))
				continue;
			int row = zScores.rowHash.get(HPOs.get(h));
			if(geneZscores == null)
			{
				geneZscores = zScores.getRow(row);
				continue;
			}
			
			for(int c = 0; c < zScores.cols(); c++)
			{
				int geneRow = geneZscores.rowHash.get(zScores.getColHeaders()[c]); 
				geneZscores.matrix.add(geneRow, 0, zScores.matrix.get(row, c));
			}	
		}
		if(geneZscores == null)
			return -1;
		geneZscores.sortCol(0);
		rank = geneZscores.rowHash.get(gene)+1;
		return rank;
	}

	private static String convertGStoEng(int col, String fileName) throws Exception 
	{
		String writeName = fileName.replace(".txt", "_ENSG.txt");
		Hashtable<String,String> geneSymbolToEnsembl = readConversionHash(); 
		BufferedReader reader = new BufferedReader(new FileReader(new File(fileName)));
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(writeName)));
		String line = reader.readLine();//skip the header line
		while((line = reader.readLine())!= null)
		{
			String[] eles = line.split("\t");
			String geneSymbol = eles[col];
			if(geneSymbolToEnsembl.get(geneSymbol) == null)//if it doesn't exist in the geneSymbolToEnsemblFN just skip it.
				continue;
			line = line.replace(geneSymbol, geneSymbolToEnsembl.get(geneSymbol));
			//System.out.println(line);
			writer.write(line+"\n");
		}
		
		writer.close();
		reader.close();
		return writeName;
	}

	private static Hashtable<String, String> readConversionHash() throws Exception 
	{
		Hashtable<String,String> geneSymbolToEnsembl = new Hashtable<String,String>();
		BufferedReader reader = new BufferedReader(new FileReader(new File(geneSymbolToEnsemblFN)));
		String line = null;
		while((line = reader.readLine()) != null)
		{
			String[] eles = line.split("\t");
			String geneSymbol = eles[1];
			String geneEnsgID = eles[0];
			geneSymbolToEnsembl.put(geneSymbol, geneEnsgID);
		}	
		reader.close();
		return geneSymbolToEnsembl;
	}

	private static Hashtable<String, ArrayList<String>> readOMIM_HPOs(String iDS_ConvertedFN) throws IOException {
		Hashtable<String, ArrayList<String>> OMIM_HPO = new Hashtable<String, ArrayList<String>>();
		BufferedReader reader = new BufferedReader(new FileReader(new File(iDS_ConvertedFN)));
		String line = reader.readLine();//skip the header line
		while((line = reader.readLine())!= null)
		{
			String[] eles = line.split("\t");
			String OMIM = eles[0];
			String HPO = eles[3];
//			System.out.println(OMIM+"\t"+ HPO);
			addToTable(OMIM,HPO, OMIM_HPO);
		}
		reader.close();
		return OMIM_HPO;
	}

	private static void addToTable(String OMIM, String HPO, Hashtable<String, ArrayList<String>> OMIM_HPO) 
	{
		ArrayList<String> OMIM_HPOs = OMIM_HPO.get(OMIM);
		if(OMIM_HPOs == null)
			OMIM_HPOs = new ArrayList<String>();
		OMIM_HPOs.add(HPO);
		OMIM_HPO.put(OMIM, OMIM_HPOs);
	}
	static void checkArgs(String[] args) 
	{
		if(args.length < 4)
		{
			System.out.println("Wrong arguments:\n"
					+ "1. zscoreFN=<zScoresGenesVSHPOFN.txt> - file with the zScores per gene per HPO\n"
					+ "2. dp=<diseaseToPhenotype.txt> - has OMIM terms on the 1st column and HPO terms on the 4th\n"
					+ "3. idToEnsg=<ENSGToGeneNameHGNCBiotypeChromosomeStartStopStrandAndDescriptionV75.txt> - IDs on the .. column and the ENSGids on 2nd\n"
					+ "4. clinvar=<ClinVar_variant_summary.txt> - has IDs on 5th column, PhenotypeIDs (OMIM) on 10th");
			System.exit(1);
		}
		
		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase())
			{
				case "zscorefn":
					zScoresGenesVSHPOFN =value;
				break;
				case "dp":
					diseaseGeneFN =value;
					break;
				case "idtoensg":
					geneSymbolToEnsemblFN =value;
					break;
				case "clinvar":
					clinvar =value;
					break;
				default:
					System.out.println("Incorrect argument supplied:\n"+ args[a] +"\nexiting");
					System.exit(1);
			}
		}
	}
}
