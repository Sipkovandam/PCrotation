package GeneNetwork;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;

import MatrixScripts.MyMatrix;
import Tools.FileUtils;
import Tools.Script;

public class WebsiteMatrixCreator extends Script<WebsiteMatrixCreator>
{
	//Takes a Patrick file (GeneXHPOterm with 0 or 1 to indicate whether the gene is part of the term)
	//and converts it to a Juha input file (per line gene\tTerm4\tTerm10) for this script ~/gene-network/data/populateGenesetDBTXT.js
	
	String termsUrlFns = null;///comma separated list of files. Filecontain term IDs (eg R-HSA-1059683) in first column, webUrl in second column and pathway name in third (e.g. Interleukin-6 signaling) // significant pathways only
	
	String geneXTermFns = null;//comma separated list of files. Files contain genes on rows,terms on columns. Values[][] are 1 if gene is part of term, 0 if not.
	
	String aucFns = null;//comma separated list of files. Files contain term names in first column, p-value in second, AUC in third column.
	
	String zScoreFns = null;//comma separated list of files. Files contain genes on rows,terms on columns. Values are z-scores
	
	boolean toUpperCase = true;
	
	transient MyMatrix sortOrder = null;
	transient HashMap<String,Integer> termName_To_RowNumber = null;

	@Override
	public void run()
	{
		try
		{
			HashMap<String,String> termId_To_TermName = createConversionHash(termsUrlFns);
			log("Conversion hash contains\t" + termId_To_TermName.size() + " IDs");
			
			//get the order all matrixes should be in
			termName_To_RowNumber = getTermNameToRowNumberHash(termId_To_TermName);
			
			convertPathWayUrlFiles(termId_To_TermName);
			
//			convertGeneXHpoTermFiles(termId_To_TermName);
//			
//			convertAucFns(termId_To_TermName);
//			
//			convertGeneNamesInZscoreFiles(termId_To_TermName);
			
			
			
		}catch(Exception e){e.printStackTrace();}
	}

	private void convertPathWayUrlFiles(HashMap<String, String> termId_To_TermName) throws IOException
	{
		log("Converting file:\t" + termsUrlFns);
		String outputFn = FileUtils.removeExtention(termsUrlFns)+ "_gnInputFormat.txt";
		
		BufferedReader urlFileReader = FileUtils.createReader(termsUrlFns);
		
		String[] writeBuffer = new String[this.termName_To_RowNumber.size()];//foutje, hier staan genen op de rownames 
		
		String inputLine=null;
		StringBuilder writeLine = new StringBuilder(); 
		
		while((inputLine=urlFileReader.readLine())!=null)
		{
			inputLine=inputLine.replace("\"", "");
			writeLine.setLength(0);
			String[] eles =inputLine.split("\t");
			String term = eles[2].toUpperCase();
			writeLine.append(term);
			writeLine.append("\t-\t");
			writeLine.append(eles[1]);
			writeLine.append("\t-");;
			
			
			if(termName_To_RowNumber.get(term)!=null)
			{
				
				int rowNumber = termName_To_RowNumber.get(term);
				writeBuffer[rowNumber]=writeLine.toString();
			}
		}
		
		writeLineBuffer(writeBuffer, outputFn);			
		
		log("File written to:\t" + outputFn);
	}

	private HashMap<String, Integer> getTermNameToRowNumberHash(HashMap<String, String> termId_To_TermName)
	{

		String[] termNames = termId_To_TermName.values().toArray(new String[termId_To_TermName.values().size()]);
		Arrays.sort(termNames);
		HashMap<String, Integer> termName_To_RowNumber = new HashMap<String, Integer>();
		for(int t = 0; t < termNames.length; t++)
			termName_To_RowNumber.put(termNames[t].toUpperCase(), t);
		
		return termName_To_RowNumber;
	}

	private void convertGeneNamesInZscoreFiles(HashMap<String, String> termId_To_TermName)
	{
		if(zScoreFns==null)
		{
			log("No zScoreFns files supplied, skipping this step");
			return;
		}
		
		String[] filesToConvert = zScoreFns.split(",");
		
		for(String fn : filesToConvert)
		{
			MyMatrix zscoreMatrix = new MyMatrix(fn);
			for(int c = 0; c < zscoreMatrix.cols(); c++)
			{
				String termName = termId_To_TermName.get(zscoreMatrix.colNames[c]);
				zscoreMatrix.colNames[c]=termName;
			}
			String writeFn = FileUtils.removeExtention(fn)+"_termNames.txt.gz";
			zscoreMatrix.transpose();
			zscoreMatrix.keepIDs(this.termName_To_RowNumber);
			zscoreMatrix.transpose();
			zscoreMatrix.write(writeFn);
			log("File written to: \t" + writeFn);
		}
	}

	private void convertAucFns(HashMap<String, String> termId_To_TermName) throws FileNotFoundException, IOException
	{
		if(aucFns==null)
		{
			log("No AUC files supplied, skipping this step");
			return;
		}
		
		String[] filesToConvert = aucFns.split(",");
		String header = "Term\tP_value\tAUC\tnr_genes";
		
		for(String fn : filesToConvert)
		{
			log("Converting the following file:\t" + fn);
			String outputFn = FileUtils.removeExtention(fn)+ "_gnInputFormat.txt";
			MyMatrix aucMatrix = new MyMatrix(fn, true, false);
			String[] writeBuffer = new String[this.termName_To_RowNumber.size()]; 

			for(int r = 0; r < aucMatrix.rows(); r++)
			{
				String termId=aucMatrix.rowNames[r];
				String termName = termId_To_TermName.get(termId);
				if(termName==null)
					continue;
				if(toUpperCase)
					termName=termName.toUpperCase();
				log("termname=\t" + termName);
				String line = termName + "\t" + aucMatrix.values[r][1] + "\t" + aucMatrix.values[r][2] + "\t"+ aucMatrix.values[r][0];
				int rowNumber = termName_To_RowNumber.get(termName);
				writeBuffer[rowNumber]=line.toString();
			}
			writeLineBuffer(writeBuffer, outputFn, header);
			log("Reformattted AUC file written to:\t" + outputFn);
		}	
	}

	private void convertGeneXHpoTermFiles(HashMap<String, String> termId_To_TermName) throws FileNotFoundException, IOException
	{
		if(geneXTermFns==null)
		{
			log("No geneXHpoTermFns files supplied, skipping this step");
			return;
		}
		
		String[] filesToConvert = geneXTermFns.split(",");
		
		for(String fn: filesToConvert)
		{
			log("Converting file:\t" + fn);
			String outputFn = FileUtils.removeExtention(fn)+ "_gnInputFormat.txt";
			
			MyMatrix geneInTerm = new MyMatrix(fn);
			
			//String[] writeBuffer = new String[this.termName_To_RowNumber.size()];//foutje, hier staan genen op de rownames 
			String[] writeBuffer = new String[geneInTerm.rows()];
			
			StringBuilder line = new StringBuilder(); 
			for(int r = 0; r < geneInTerm.rows(); r++)
			{
				line.setLength(0);
				line.append(geneInTerm.rowNames[r]);
				
				for(int c = 0; c < geneInTerm.cols(); c++)
				{
					if(geneInTerm.values[r][c]==1)
					{
						
						String termName=termId_To_TermName.get(geneInTerm.colNames[c]);
						if(termName==null)
							continue;
						
						line.append("\t");
						if(toUpperCase)
							termName=termName.toUpperCase();
						line.append(termName);
					}
				}
				log("geneTerm=\t" + geneInTerm.rowNames[r]);
				//int rowNumber = termName_To_RowNumber.get(geneInTerm.rowNames[r]);
				writeBuffer[r]=line.toString();
			}
			
			writeLineBuffer(writeBuffer, outputFn);			
			
			log("File written to:\t" + outputFn);
		}
	}
	private void writeLineBuffer(String[] writeBuffer, String outputFn) throws FileNotFoundException, IOException
	{
		writeLineBuffer(writeBuffer, outputFn,null);	
	}
	
	private void writeLineBuffer(String[] writeBuffer, String outputFn, String header) throws FileNotFoundException, IOException
	{
		BufferedWriter geneToTermsWriter = FileUtils.createWriter(outputFn);
		
		if(header!=null)
			geneToTermsWriter.write(header+"\n");

		for(String line : writeBuffer)
			geneToTermsWriter.write(line.toString()+"\n");

		geneToTermsWriter.close();
	}

	private HashMap<String, String> createConversionHash(String termIdToTermNameFns)
	{
		HashMap<String, String> termId_To_TermName = new HashMap<String, String>();
		String[] fns = termIdToTermNameFns.split(",");
		for(String fn : fns)
		{
			termId_To_TermName.putAll(FileUtils.readStringStringHash(fn, 0, 2,"noQuotes"));
		}
		for(String termId : termId_To_TermName.keySet())
		{
			String termName = termId_To_TermName.get(termId);
			termId_To_TermName.put(termId, termName.toUpperCase());	
		}
		
		return termId_To_TermName;
	}

}
