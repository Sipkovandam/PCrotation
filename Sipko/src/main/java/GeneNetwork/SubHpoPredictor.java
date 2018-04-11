package GeneNetwork;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import MatrixScripts.MyMatrix;
import Tools.FileUtils;
import Tools.Script;

public class SubHpoPredictor extends Script<SubHpoPredictor>
{
	final String webUrl="https://mseqdr.org/hpo_browser.php?";
	String ensemblToGeneSymbolFn="E:/Groningen/Data/Annotation/GRCh37/ENSGToGeneNameV75.txt";
		
	String geneHpoZscoreMatrixFn = "E:/Groningen/Test/GeneNetwork.SubHpoPredictor/HpoZscores/hpo_predictions.txt.gz";	
	String dnaNumberToCausalGeneFn = "E:/Groningen/Test/GeneNetwork.SubHpoPredictor/dnaNumberToLocation.txt";
	String writeFolder = "E:/Groningen/Test/GeneNetwork.SubHpoPredictor/results/";
	boolean includeSubHpoTermAnalysis=true;
	
	transient HashMap<String,String> ensembl_To_GeneSymbol = FileUtils.readStringStringHash(ensemblToGeneSymbolFn);
	transient HashMap<String,String> geneSymbol_To_ensembl = FileUtils.readStringStringHash(ensemblToGeneSymbolFn,1,0);
	
	@Override
	public void run()
	{
		try
		{
			writeFolder=FileUtils.makeFolderNameEndWithSlash(writeFolder);
			FileUtils.makeDir(writeFolder);
			
			HashMap<String,String> dnaNumber_To_HpoTerm = FileUtils.readStringStringHash(dnaNumberToCausalGeneFn,0,1);
			HashMap<String,String> dnaNumber_To_CausalGenes = FileUtils.readStringStringHash(dnaNumberToCausalGeneFn,0,2);
			HashMap<String,String> dnaNumber_To_FileLocation = FileUtils.readStringStringHash(dnaNumberToCausalGeneFn,0,3);

			BufferedWriter faultyGeneSymbolsWriter = FileUtils.createWriter(writeFolder+"incorrectGeneSymbols.txt");
			
			BufferedWriter causalGeneRankWriter = FileUtils.createWriter(writeFolder+"summary.txt");
			String header = "DNAnumber\tcausalGene\trank\ttotal\textraInfo\n";
			causalGeneRankWriter.write(header);

			MyMatrix geneHpoZscoreMatrix = new MyMatrix(geneHpoZscoreMatrixFn);

			for(String dnaNumber : dnaNumber_To_CausalGenes.keySet())
			{
				//if any info is missing skip this DNA number
				if(!FileUtils.hasValueForEntry(dnaNumber_To_HpoTerm, dnaNumber) || !FileUtils.hasValueForEntry(dnaNumber_To_CausalGenes, dnaNumber) || !FileUtils.hasValueForEntry(dnaNumber_To_FileLocation, dnaNumber))
					continue;
				
				String[] hpoTerms = dnaNumber_To_HpoTerm.get(dnaNumber).replace("\"", "").split(",");
				
				String gavinPathogenicFilteredFn=dnaNumber_To_FileLocation.get(dnaNumber);
				
				//if GAVIN pathogenic file does not exist skip to the next dnaNumber
				if(!new File(gavinPathogenicFilteredFn).exists())
					continue;
				
				String writeFn = writeFolder+dnaNumber+"_hpoRanks.txt";
				String[] causalGenes = dnaNumber_To_CausalGenes.get(dnaNumber).split("/");	
				
				log("Running sample:\t "+dnaNumber);
				
				runHpoSubtermAnalysis(gavinPathogenicFilteredFn, hpoTerms, writeFn, dnaNumber, causalGenes,causalGeneRankWriter, geneHpoZscoreMatrix,faultyGeneSymbolsWriter);
			}
			faultyGeneSymbolsWriter.close();
			causalGeneRankWriter.close();
		}catch(Exception e){e.printStackTrace();}
	}

	private void runHpoSubtermAnalysis(String gavinPathogenicFilteredFn, String[] hpoTerms, String writeFn, String dnaNumber, String[] causalGenes, BufferedWriter causalGeneRankWriter, MyMatrix geneHpoZscoreMatrix, BufferedWriter faultyGeneSymbolsWriter) throws IOException
	{
		
		//get parent HPO z-scores
		
		//header line
		String[] header =new String[]{"GeneName"};
		String[] outputLines=null;

		double[] parentZscoreSum = new double[geneHpoZscoreMatrix.rows()];
		int nParents = 0;
		double[] weightedParentNewZscoreSum=new double[geneHpoZscoreMatrix.rows()];
		
		for(String hpoTerm : hpoTerms)
		{
			if(geneHpoZscoreMatrix.getColValues(hpoTerm)==null)
			{
				log("HPOterm:\t" + hpoTerm + "\tnot present in HPOzscoreMatrix\t" + geneHpoZscoreMatrixFn);
				continue;
			}
			
			//get z-score for parent for each gene
			double[] hpoZ_scores = geneHpoZscoreMatrix.getColValues(hpoTerm);
			header[0]+="\tParent_"+hpoTerm;
			if(outputLines==null)
				outputLines=intitiateLines(hpoZ_scores, geneHpoZscoreMatrix);
			else
			{
				addZscoresToLines(outputLines,hpoZ_scores);
			}
			
			parentZscoreSum=FileUtils.sumDoubleArrays(hpoZ_scores, parentZscoreSum);
			
			if(includeSubHpoTermAnalysis)
				weightedParentNewZscoreSum=childHpoAnalysis(hpoTerm, geneHpoZscoreMatrix, header, outputLines,weightedParentNewZscoreSum, hpoZ_scores);
			nParents++;
		}
		
		if(outputLines==null)
		{
			log("No HPO terms or child HPO terms found skipping:\t" + dnaNumber);
			return;
		}
		
		//calculate weighted z-score for each gene for the original parent terms
		parentZscoreSum=weightZscores(parentZscoreSum,nParents);
		header[0]+="\tweightedParents";
		addZscoresToLines(outputLines,parentZscoreSum);
		
		if(includeSubHpoTermAnalysis)
		{
			//calculate weighted z-score for each gene for the new weighted parent scores
			weightedParentNewZscoreSum=weightZscores(weightedParentNewZscoreSum,nParents);
			header[0]+="\tnewZscoresWeighted_Parents";
			addZscoresToLines(outputLines,weightedParentNewZscoreSum);
		}
		
		BufferedWriter scoreWriter = FileUtils.createWriter(writeFn);
		scoreWriter.write(header[0]+"\n");
		for(String line: outputLines)
			scoreWriter.write(line+"\n");
		
		scoreWriter.close();

		log("All genes HPO prediction file written to:\t" + writeFn);
		
		//write only genes containing pathogenic variants
		String writeFnPathogenicOnly = FileUtils.removeExtention(writeFn)+"_pathogenicOnly.txt";
		BufferedWriter pathoGenicScoreWriter = FileUtils.createWriter(writeFnPathogenicOnly);
		pathoGenicScoreWriter.write(header[0]+"\n");
		
		HashSet<String> pathogenicGenes = FileUtils.readHashSet(gavinPathogenicFilteredFn,8);
		for(String line: outputLines)
		{
			String ensemblName = line.split("\t")[0];//ontzettend lelijk, ik weet het
			String geneSymbol = ensembl_To_GeneSymbol.get(ensemblName);
			if(pathogenicGenes.contains(geneSymbol))
				pathoGenicScoreWriter.write(line+"\n");
			
		}
		pathoGenicScoreWriter.close();
	
		//write the rank that the causal gene has in the genenetwork+GAVIN results, to a separate file.
		MyMatrix file = new MyMatrix(writeFnPathogenicOnly);
		file.sortCol(file.cols()-1);
		file.write(writeFnPathogenicOnly);		
		
		int causalGeneRow = -2;
		for(String causalGene:causalGenes)
		{
			String ensemblId=geneSymbol_To_ensembl.get(causalGene);
			if(ensemblId==null)
			{
				log("Warning no ensemblID found for: \t" + causalGene);
				faultyGeneSymbolsWriter.write(causalGene+"\n");
				continue;
			}
			
			log("ensemblID=" + ensemblId + "\tdnaNumber=\t" +dnaNumber);
			if(ensemblId!= null && file.getRowHash().containsKey(ensemblId))
				causalGeneRow=file.getRowHash().get(ensemblId);
			String summaryLine = dnaNumber+"\t"+causalGene+"\t"+(causalGeneRow+1)+"\t"+file.rows();
			if(causalGeneRow>=0)
				summaryLine+="\t"+file.getRowInfoInOneLine(causalGeneRow);
			
			causalGeneRankWriter.write(summaryLine +"\n");
		}
		
		log("Genes with pathogenic mutations only file written to:\t" + writeFnPathogenicOnly);

	}

	private double[] childHpoAnalysis(String hpoTerm, MyMatrix geneHpoZscoreMatrix, String[] header, String[] outputLines, double[] weightedParentNewZscoreSum, double[] parentZscore) throws IOException
	{
		String[] hpoChildren = getChildren(hpoTerm);
		
		//get z-score for each child for each gene
		double[] zScoreSum = new double[geneHpoZscoreMatrix.rows()]; 
		int n = 0;
		for(String hpoChild : hpoChildren)
		{
			log("getting z-score for HPO term:\t" + hpoChild);
			if(geneHpoZscoreMatrix.getColHash().get(hpoChild)== null)
				continue;
			
			int hpoColumn = geneHpoZscoreMatrix.getColHash().get(hpoChild);
			double[] hpoZ_scoresChild = geneHpoZscoreMatrix.getColValues(hpoColumn);
			
			//add to output file
			header[0]+="\tChild_"+hpoChild;
			addZscoresToLines(outputLines,hpoZ_scoresChild);
			
			zScoreSum=FileUtils.sumDoubleArrays(hpoZ_scoresChild, zScoreSum);
			n++;
		} 
		
		header[0]+="\tweighted_"+hpoTerm;
		if(n>0)
		{
			//calculate weighted z-score for each gene
			zScoreSum=weightZscores(zScoreSum,n);
		}
		else//non of the HPO child terms where significant and thus non where included
		{
			zScoreSum=parentZscore;
		}
		
		addZscoresToLines(outputLines,zScoreSum);
		
		weightedParentNewZscoreSum=FileUtils.sumDoubleArrays(weightedParentNewZscoreSum, zScoreSum);
		return weightedParentNewZscoreSum;
	}

	private double[] weightZscores(double[] zScoreSum, int n)
	{
		double divider = Math.pow(n, 0.5);
		for(int r = 0; r < zScoreSum.length; r++)
		{
			zScoreSum[r]/=divider;
		}
		return zScoreSum;
	}

	private void addZscoresToLines(	String[] outPutlines,
									double[] hpoZ_scoresChild)
	{
		for(int r = 0; r < hpoZ_scoresChild.length;r++)
		{
			outPutlines[r]+="\t"+hpoZ_scoresChild[r];
		}
		
	}

	private String[] intitiateLines(double[] hpoZ_scores,
									MyMatrix geneHpoZscoreMatrix)
	{
		String[] outPutlines = new String[hpoZ_scores.length];
		
		for(int r = 0; r < hpoZ_scores.length; r++)
		{
			outPutlines[r]=geneHpoZscoreMatrix.rowNames[r]+"\t"+hpoZ_scores[r];
		}
		return outPutlines;			
	}

	private String[] getChildren(String hpoTerm) throws IOException
	{
		String hpoShort = hpoTerm.replaceAll("HP:[0]{0,8}", "");
		
		String hpoUrl = webUrl+hpoShort+";";
		log("hpoUrl=\t" + hpoUrl);
		BufferedReader webReader = FileUtils.getWebsiteReader(hpoUrl);
		
		ArrayList<String> childTerms = new ArrayList<String>();
		String line=null;
		
		//TO DO: get child terms...
		while((line=webReader.readLine())!=null)
		{
			if(!line.contains("Child Nodes:"))
				continue;
			
			String[] childLines = line.split("\\.\\.\\.\\.\\.\\.\\.\\.<a href");
			
			for(String childLine:childLines)
			{
				if(!childLine.contains("term?id="))
					continue;
				String childHpo = getHpoName(childLine);
				log("childHpo=\t" + childHpo);
				childTerms.add(childHpo);
			}	
		}
		
		return (String[]) childTerms.toArray(new String[childTerms.size()]);
	}

	private String getHpoName(String childLine)
	{
		String childHpo = childLine.split("term\\?id=")[1].split("\\\" onClick=")[0];
		return childHpo;
	}
	
}
