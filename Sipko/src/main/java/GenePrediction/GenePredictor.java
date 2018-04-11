package GenePrediction;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import MatrixScripts.MyMatrix;
import Tools.FileUtils;
import Tools.Script;

public class GenePredictor extends Script<GenePredictor>
{
	///
	
	//input variables
	String methodComment = "method can be one of the following options: getDiseaseZscores_Standard, getCausalGeneRanks";
	String method = "getDiseaseZscores_Standard";
	String orphanetFn = "E:/Groningen/Students/Patrick/Annotation/ALL_SOURCES_FREQUENT_FEATURES_diseases_to_genes_to_phenotypes.txt";
	String zScoreMatrixFn = "E:/Groningen/Test/GenePrediction.GenePredictor/hpo_predictions_bonSigOnly.txt.gz";
	
	String writeFn = "E:/Groningen/Test/GenePrediction.GenePredictor/zScoreGenesPerDisease.txt";
	String writeFnRankMatrix = "E:/Groningen/Test/GenePrediction.GenePredictor/zScoreGenesPerDiseaseRanks.txt";
	String writeFnDiseaseGenePredictions = "E:/Groningen/Test/GenePrediction.GenePredictor/diseaseGenePredictions.txt";

	//globals
	transient MyMatrix zScoreMatrix = null;
	transient OrpahNetInfo orpahNetInfo = null;
	
	@Override
	public void run()
	{
		try
		{
			log("Method =" +  method);
			switch(method)
			{
				case "getDiseaseZscores_Standard"://get matrix of genesXdiseases with the z-score of each gene for that ophanet disease
					createDiseaseZscoreMatrix();
				case "getCausalGeneRanks"://get the rank of the causal gene for each disease.
					getGeneRanks();
			}
		
		}catch(Exception e){e.printStackTrace();}
	}

	private void getGeneRanks() throws FileNotFoundException, IOException
	{
		if(new File(writeFnRankMatrix).exists())
			this.createDiseaseZscoreMatrix();
		
		MyMatrix ranksPerDisease = new MyMatrix(writeFnRankMatrix);
		
		Set<String> diseases = orpahNetInfo.disease_To_Genes.keySet();
		int c = 0;
		
		BufferedWriter diseaseGeneRankWriter= FileUtils.createWriter(writeFnDiseaseGenePredictions);
		for(String disease : diseases)
		{
			log("Disease =\t" + disease);
			//TODO: Pick parent term if not found
			if(orpahNetInfo.disease_To_HpoTerms.containsKey(disease)) //if this hpo term is not significant or does not exist
				continue;
			
			switch(method)
			{
				case "getCausalGeneRanks":
					HashSet<String> genes = orpahNetInfo.disease_To_Genes.get(disease);
					writeGeneRanks(diseaseGeneRankWriter, genes, disease, ranksPerDisease);
					break;
			}
			

			
		}
		diseaseGeneRankWriter.close();
	}

	private void writeGeneRanks(BufferedWriter diseaseGeneRankWriter,
								HashSet<String> genes,
								String disease, MyMatrix ranksPerDisease) throws IOException
	{
		int column = ranksPerDisease.getColHash().get(disease);//I dont think this can be null at this point, but might need expection if I am wrong
		for(String gene: genes)
		{
			
			int row = ranksPerDisease.getRowHash().get(gene);
			String rank = Integer.toString((int)ranksPerDisease.values[row][column]);
			
			String line = disease+"\t"+gene+"\t"+rank+"\n";
			diseaseGeneRankWriter.write(line);
		}
	}

	private void createDiseaseZscoreMatrix() throws FileNotFoundException, IOException
	{
		init();
		
		log("Input read, z-score calculation");
		MyMatrix summedZscores = new MyMatrix(zScoreMatrix.rows() ,2);	
		summedZscores.rowNames=zScoreMatrix.rowNames;

		Set<String> diseases = orpahNetInfo.disease_To_Genes.keySet();
		int c = 0;
		for(String disease : diseases)
		{
			log("Disease =\t" + disease);
			//TODO: Pick parent term if not found
			if(orpahNetInfo.disease_To_HpoTerms.containsKey(disease)) //if this hpo term is not significant or does not exist
				continue;
			
			summedZscores.setColHeader(c, disease);
			summedZscores = getZscoreMatrix(summedZscores, disease, c);
			c++;
			break;
		}
	
		log("DiseaseZscore file written to:\t " + writeFn);
		summedZscores.write(writeFn);
		
		log("Creating ranking matrix");
		//sort
		MyMatrix ranksMatrix = convertToRankMatrix(summedZscores);
		
		log("Ranks file written to:\t " + writeFnRankMatrix);
		ranksMatrix.write(writeFnRankMatrix);
	}

	private MyMatrix convertToRankMatrix(MyMatrix summedZscores)
	{
		MyMatrix ranksMatrix = new MyMatrix(zScoreMatrix.rows() ,orpahNetInfo.disease_To_Genes.size());
		ranksMatrix.rowNames= FileUtils.makeDeepCopy(zScoreMatrix.rowNames);

		ranksMatrix.colNames=summedZscores.colNames;
		
		for(int c =0; c < summedZscores.cols(); c++)
		{
			log("Col=" + c + "/" + summedZscores.cols());
			summedZscores.sortCol(c);
			
			for(int r = 0; r < ranksMatrix.rows();r++)
			{
				ranksMatrix.values[r][c]= summedZscores.getRowHash().get(ranksMatrix.rowNames[r]);
			}
		}
		
		return ranksMatrix;
	}

	private void init() throws FileNotFoundException, IOException
	{
		orpahNetInfo=new OrpahNetInfo(orphanetFn);
		zScoreMatrix=new MyMatrix(zScoreMatrixFn);
		writeFnRankMatrix=FileUtils.removeExtention(writeFn)+"_ranks.txt";
	}

	private MyMatrix getZscoreMatrix(MyMatrix summedZscores,
									String disease, int c)
	{
		HashSet<String> hpoTerms = orpahNetInfo.disease_To_HpoTerms.get(disease);

		switch(method)
		{
			case "getDiseaseZscores_Standard:":
				summedZscores = getWeightedZscores(hpoTerms, summedZscores, c);
				
		}
		return summedZscores;
	}

	private MyMatrix getWeightedZscores(HashSet<String> hpoTerms, MyMatrix summedZscores, int c)
	{
		double nHpoTerms = sumZscores(hpoTerms, summedZscores, c);
		
		double denominator = java.lang.Math.pow(nHpoTerms, 0.5);
		
		for(int r = 0; r < zScoreMatrix.rows();r++)
		{
			summedZscores.values[r][c]/=denominator;
		}
		
		return summedZscores;
	}

	private double sumZscores(HashSet<String>  hpoTerms, MyMatrix summedZscores, int c)
	{
		double nHpoTerms = 0;
		for(String hpoTerm : hpoTerms)
		{
			if(!zScoreMatrix.getColHash().contains(hpoTerm))
				log("Warning this HPO term is not present in the z-score matrix:\t" + hpoTerm);
			
			for(int r = 0; r < zScoreMatrix.rows();r++)
			{
				int zscoreCol = zScoreMatrix.getColHash().get(hpoTerm);
				
				summedZscores.values[r][c]+=zScoreMatrix.values[r][zscoreCol];
			}
			nHpoTerms++;
		}
		return nHpoTerms;
	}
}
