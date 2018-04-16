package GeneNetwork;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import MatrixScripts.MyMatrix;
import Tools.FileUtils;
import Tools.Script;

public class WebsiteMatrixCreator extends Script<WebsiteMatrixCreator>
{
	//Takes a Patrick file (GeneXHPOterm with 0 or 1 to indicate whether the gene is part of the term)
	//and converts it to a Juha input file (per line gene\tTerm4\tTerm10) for this script ~/gene-network/data/populateGenesetDBTXT.js

	String termsUrlFn = null;///comma separated list of files. File contain term IDs (eg R-HSA-1059683) in first column, webUrl in second column and pathway name in third (e.g. Interleukin-6 signaling) // significant pathways only

	String geneXTermFn = null;//comma separated list of files. Files contain genes on rows,terms on columns. Values[][] are 1 if gene is part of term, 0 if not.

	String aucFn = null;//comma separated list of files. Files contain term names in first column, p-value in second, AUC in third column.

	String zScoreFn = null;//comma separated list of files. Files contain genes on rows,terms on columns. Values are z-scores

	String writeFolder = null;
	
	boolean toUpperCase = true;

	transient MyMatrix sortOrder = null;
	transient HashMap<String, Integer> term_To_RowNumber = null;

	@Override
	public void run()
	{
		try
		{
			FileUtils.makeDir(writeFolder);
			
			HashMap<String, String> termId_To_TermName = createConversionHash(termsUrlFn);
			log("Conversion hash contains\t" + termId_To_TermName.size() + " IDs");

			//get the order all matrixes should be in
			term_To_RowNumber = getTermToRowNumberHash(termId_To_TermName);

			convertPathWayUrlFiles(termId_To_TermName);

			convertGeneXHpoTermFiles(termId_To_TermName);
			//			
			convertAucFns(termId_To_TermName);
			//			
			convertGeneNamesInZscoreFiles(termId_To_TermName);

		} catch (Exception e)
		{
			e.printStackTrace();
		}
	}

	private void convertPathWayUrlFiles(HashMap<String, String> termId_To_TermName) throws IOException
	{
		log("Converting file:\t" + termsUrlFn);
		String outputFn = FileUtils.makeFolderNameEndWithSlash(writeFolder)+new File(FileUtils.removeExtention(termsUrlFn) + "_gnInputFormat.txt").getName();

		BufferedReader urlFileReader = FileUtils.createReader(termsUrlFn);

		int aantal = Collections.max(this.term_To_RowNumber.values()) + 1;

		String[] writeBuffer = new String[aantal];

		String inputLine = null;
		StringBuilder writeLine = new StringBuilder();

		while ((inputLine = urlFileReader.readLine()) != null)
		{
			inputLine = inputLine.replace(	"\"",
											"");
			writeLine.setLength(0);
			String[] eles = inputLine.split("\t");
			String termId = eles[0].toUpperCase();
			String termName = eles[2].toUpperCase();
			writeLine.append(termName);
			writeLine.append("\t-\t");
			writeLine.append(eles[1]);
			writeLine.append("\t-");

			if (term_To_RowNumber.get(termName) != null)
			{
				int rowNumber = term_To_RowNumber.get(termName);
				writeBuffer[rowNumber] = writeLine.toString();
			}
			else
				log("term=\t" + termName + "\tsizeofHash=\t" + term_To_RowNumber.size() + "\tfirstKey=\t" + term_To_RowNumber.keySet().iterator().next());

		}

		writeLineBuffer(writeBuffer,
						outputFn);

		log("File written to:\t" + outputFn);
	}

	private HashMap<String, Integer> getTermToRowNumberHash(HashMap<String, String> termId_To_TermName)
	{

		String[] termNames = termId_To_TermName.values().toArray(new String[termId_To_TermName.values().size()]);
		String[] termCodes = termId_To_TermName.keySet().toArray(new String[termId_To_TermName.values().size()]);
		HashMap<String, String> termName_To_termId=FileUtils.invertHash_KeysToValues_ValuesToKeys(termId_To_TermName);
		Arrays.sort(termNames);
		HashMap<String, Integer> termName_To_RowNumber = new HashMap<String, Integer>();
		for (int t = 0; t < termNames.length; t++)
		{
			termName_To_RowNumber.put(	termNames[t],
										t);
			termName_To_RowNumber.put(	termName_To_termId.get(termNames[t]),
										t);
		}
		return termName_To_RowNumber;
	}

	private void convertGeneNamesInZscoreFiles(HashMap<String, String> termId_To_TermName)
	{
		if (zScoreFn == null)
		{
			log("No zScoreFn files supplied, skipping this step");
			return;
		}

		MyMatrix zscoreMatrix = new MyMatrix(zScoreFn);
		for (int c = 0; c < zscoreMatrix.cols(); c++)
		{
			String termName = termId_To_TermName.get(zscoreMatrix.colNames[c]);
			zscoreMatrix.colNames[c] = termName;
		}
		String writeFn = FileUtils.removeExtention(zScoreFn) + "_termNames.txt.gz";
		zscoreMatrix.transpose();
		zscoreMatrix.keepIDs(this.term_To_RowNumber);
		zscoreMatrix.transpose();
		zscoreMatrix.write(writeFn);
		log("Zscore matrix file has:\t"+zscoreMatrix.cols()+" coluns, File written to: \t" + writeFn);
	}

	private void convertAucFns(HashMap<String, String> termId_To_TermName) throws FileNotFoundException, IOException
	{
		if (aucFn == null)
		{
			log("No AUC files supplied, skipping this step");
			return;
		}

		String header = "Term\tP_value\tAUC\tnr_genes";

		log("Converting the following file:\t" + aucFn);
		String outputFn = FileUtils.makeFolderNameEndWithSlash(writeFolder)+new File(FileUtils.removeExtention(aucFn) + "_gnInputFormat.txt").getName();;//FileUtils.removeExtention(aucFn) + "_gnInputFormat.txt";
		MyMatrix aucMatrix = new MyMatrix(	aucFn,
											true,
											false);

		int aantal = Collections.max(this.term_To_RowNumber.values())+1;
		String[] writeBuffer = new String[aantal];

		for (int r = 0; r < aucMatrix.rows(); r++)
		{
			String termId = aucMatrix.rowNames[r];
			String termName = termId_To_TermName.get(termId);

			if (termName == null)
				continue;
			if (toUpperCase)
				termName = termName.toUpperCase();
			String line = termName + "\t" + aucMatrix.values[r][1] + "\t" + aucMatrix.values[r][2] + "\t" + aucMatrix.values[r][0];
			int rowNumber = term_To_RowNumber.get(termName);
			writeBuffer[rowNumber] = line.toString();
		}
		writeLineBuffer(writeBuffer,
						outputFn,
						header);
		log("Reformattted AUC file written to:\t" + outputFn);

	}

	private void convertGeneXHpoTermFiles(HashMap<String, String> termId_To_TermName) throws FileNotFoundException, IOException
	{
		if (geneXTermFn == null)
		{
			log("No geneXHpoTermFns files supplied, skipping this step");
			return;
		}

		log("Converting file:\t" + geneXTermFn);
		String outputFn = FileUtils.makeFolderNameEndWithSlash(writeFolder)+new File(FileUtils.removeExtention(geneXTermFn) + "_gnInputFormat.txt").getName();//FileUtils.removeExtention(geneXTermFn) + "_gnInputFormat.txt";

		MyMatrix geneInTerm = new MyMatrix(geneXTermFn);

		//String[] writeBuffer = new String[this.termName_To_RowNumber.size()];//foutje, hier staan genen op de rownames 
		String[] writeBuffer = new String[geneInTerm.rows()];

		StringBuilder line = new StringBuilder();
		for (int r = 0; r < geneInTerm.rows(); r++)
		{
			line.setLength(0);
			line.append(geneInTerm.rowNames[r]);
			;

			for (int c = 0; c < geneInTerm.cols(); c++)
			{
				if (geneInTerm.values[r][c] == 1)
				{
					String termName = termId_To_TermName.get(geneInTerm.colNames[c]);

					if (termName == null)
						continue;

					line.append("\t");
					if (toUpperCase)
						termName = termName.toUpperCase();
					line.append(termName);
				}
			}
			//int rowNumber = termName_To_RowNumber.get(geneInTerm.rowNames[r]);
			writeBuffer[r] = line.toString();
		}

		writeLineBuffer(writeBuffer,
						outputFn);

		log("File written to:\t" + outputFn);

	}

	private void writeLineBuffer(	String[] writeBuffer,
									String outputFn) throws FileNotFoundException, IOException
	{
		writeLineBuffer(writeBuffer,
						outputFn,
						null);
	}

	private void writeLineBuffer(	String[] writeBuffer,
									String outputFn,
									String header) throws FileNotFoundException, IOException
	{
		BufferedWriter geneToTermsWriter = FileUtils.createWriter(outputFn);

		if (header != null)
			geneToTermsWriter.write(header + "\n");

		for (String line : writeBuffer)
		{
			if (line == null)
			{
				log("Warning line is:\t" + line + " in writeLineBuffer; skippng line");
				continue;
			}
			geneToTermsWriter.write(line + "\n");
		}

		geneToTermsWriter.close();
	}

	private HashMap<String, String> createConversionHash(String termIdToTermNameFn)
	{
		HashMap<String, String> termId_To_TermName = new HashMap<String, String>();

		termId_To_TermName.putAll(FileUtils.readStringStringHash(	termIdToTermNameFn,
																	0,
																	2,
																	"noBullshit"));
		
		HashMap<String, String> termId_To_TermName_upperCase = new HashMap<String, String>();
		HashSet<String> checkDoubleValues = new HashSet<String>();
		for (String termId : termId_To_TermName.keySet())
		{
			String termName = termId_To_TermName.get(termId);
			if(checkDoubleValues.contains(termName))//this is necessary because differen matrixes use different ID types and some are ambiguous meaning they lack a row/column compared to the other files
				continue;

			checkDoubleValues.add(termName);
			termId_To_TermName_upperCase.put(	termId,
									termName.toUpperCase());
		}

		return termId_To_TermName_upperCase;
	}

}
