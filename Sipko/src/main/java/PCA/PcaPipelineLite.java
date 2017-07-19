package PCA;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;

import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.concurrent.ConcurrentCorrelation;
import umcg.genetica.math.stats.concurrent.ConcurrentCovariation;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.file.Paths;
import java.text.DateFormat;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.Hashtable;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerConfigurationException;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.inference.TTest;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Text;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import JuhaPCA.PCA;
import MatrixScripts.LogTransform;
import MatrixScripts.MatrixStruct;
import MatrixScripts.MyMatrix;
import Tools.ExecCommand;
import Tools.FileUtils;
import Tools.Script;
import eqtlmappingpipeline.normalization.Normalizer;
import no.uib.cipr.matrix.NotConvergedException;

public class PcaPipelineLite extends Script<PcaPipelineLite>
{
	/**
	 * 1. load expression matrix 
	 * 2. removeDuplicates (any samples with correlation above <removeDuplicates> threshold) 
	 * 3. put genes on rows if they are not already 
	 * 4. Sort genes by chromosome and remove any that are not on in the <chromLocationsFile> 
	 * 5. add random values (add a random value of <randomAddValue> to any gene with an expression lower than <randomAddValue>) 
	 * 6. remove low expression (remove genes that do not have an expression of at least <minExpression> in at least <minSamplesExpressed> samples. 
	 * 7. remove a single gene with name <gene> (if you want to) 
	 * 8. quantile normalize 
	 * 8. Voom normalize 
	 * 8. Rlog normalization (without logging data). This is referring to part of the normalization DEseq uses 
	 * 9. 2log transform (adds <value> before doing transformation) 
	 * 10. Take the <highestExpressed> % genes 
	 * 12. keep genes with highest <TopVariance> % variance 
	 * 13. correct Input for standard deviation if <correctInputForSTdevs> == true 
	 * 14. calculate spearman correlation if <spearman> >=0. Uses <spearman> as minimum expression (if lower, expression is set to this number) 
	 * 15. if <setLowestToAverage> == True, all the lowest values to the average (if allowing log(0), all 0 will become this average)) 
	 * 16. Calculate column average and store them 
	 * 17. Correct for column averages 
	 * 18. Calculate rowaverages and store them 
	 * 19. Correct for row averages 
	 * 20. correctInputForSTdevsAfterCenter 
	 * 21. Put genes on rows if they are not already 
	 * 22. Remove genes that have no variance 
	 * 23. calculate correlation/covariance matrix 
	 * 24. PCA over matrix 
	 * 25. Calculate cronbach alphas (and therefore PCscores over samples) Different script(s): 
	 * 26. Take all PCs>0.7 and recalculate gene eigenvectors 
	 * 27. Put in database
	 **/
	// static CreateGeneEigenvectorFileVariables Var = new CreateGeneEigenvectorFileVariables();
	private static final long serialVersionUID = -8656788330691209748L;
	private String comments = "//all variables in this file containing (Comment) are comments on the variable below and do not need to be initiated";
	String runModeComment = "//0 for PCA+data correction; 1 for data correction alone";
	int runMode = 0; // 0 = PCA, 1 = datacorrection
	String expFileComment = "MANDATORY; e.g.: /root/diretory/expression.txt; //expression file on which the PCA should be conducted";
	private String expFile = "/root/diretory/expression.txt"; // expression file

	String writeFolderComment = "OPTIONAL; e.g.: /root/diretory/pcaResults/; //folder where all output should be written (default: removeExtention(expFile)+/)";
	private String writeFolder = null; // Folder where to write

	String genesToIncludeComment = "OPTIONAL; e.g.: /root/diretory/includeGenes.txt; //A file containing a list of genes you wish to include. Any genes not in this file are excluded. Should have gene names in a column (first column). First row is excluded (header)";
	private String genesToInclude = null; // List of genes. This list is selected prior to any other steps. Other genes are discarded

	String isPcaOverGenesComment = "//true for PCA over genes, false for PCA over samples";
	private boolean isPcaOverGenes = true; // if false it does the PCA over the samples rather than the genes
	String correlationComment = "//true for correlation, false for covariance";
	private boolean correlation = false; // if false uses covariance

	String centerGenesComment = "//true centers the expression for the genes (sets average to 0); should be used when conducting PCA over genes. If PCA over genes, use centerSamples=false";
	private boolean centerGenes = true; // adjusts the expression of genes so that the average expression of a sample become 0 (centering over the samples)
	String centerSamplesComment = "//true centers expression over the samples (sets average to 0); should be used when conducting PCA over samples. Centering over samples happens after centering over genes. If PCA over samples, reccommended to also use centerGenes=true";
	private boolean centerSamples = false; // adjusts the expression of genes so that the average expression of a sample become 0 (centering over the samples)

	String log2Comment = "//true for logging the data before conducting the PCA";
	private boolean log2 = true;
	String addLogValComment = "//Value to add before taking the logaritm when logging data (if <log2>==true) before conducting the PCA. Should be larger than 0";
	private double addLogVal = 0.5; // Value to add before taking the logaritm
	String deSeqNormComment = "//true for DESeq normalization of the data before the PCA";
	private boolean deSeqNorm = true; // If true uses DEseq normalization

	// if these variables are set the PC correction will be run after the eigenvectors are created (or if eigenvectors already exist only this part is run)
	String sampleFileComment = "Optional; e.g. /root/directory/filesToCorrect.txt//Filename of the file containing the samples for which the expression should be corrected";
	private String sampleFile = null;
	String PCsComment = "Optional; e.g. 1,4,7-10  //Princial components to correct for function takes format like: 1,4,5-10 or 3-10";
	private String PCs = null;// function takes format like: "1,4,5-10,3-10"
	String writeFolderCorrectedComment = "Optional; e.g. /root/directory/corrected/  //folder where the files corrected for principal components should be written (default=removeExtention(sample.getName())+/";
	private String writeFolderCorrected = "";
	String avgStdevFolderComment = "Optional; e.g. /root/directory/PCAresults/  //this is the folder with the files for the z-score calculations containing the averages and standard deviations to be used for this. allows z-scores to be calculated using the BBMRI averages and stdevs";
	private String avgStdevFolder = null; // this is the folder with the files for the z-score calculations containing the averages and standard deviations to be used for this. allows z-scores to be calculated using the BBMRI averages and stdevs
	private transient double zScoresCutoff = Double.parseDouble("0");
	private transient boolean correctResultsForSTdevs = false;
	private transient String tempName;

	public void run()
	{
		try
		{
			this.run(this.runMode);
		} catch (Exception e)
		{
			e.printStackTrace();
		}
	}

	public void run(int runMode) throws IOException, NotConvergedException, InterruptedException, ParserConfigurationException, TransformerException
	{
		this.writeConfig(	this.jsonFN,
							this,
							false,
							true);
		switch (runMode)
		{
		case 0:
		{
			if (this.writeFolder == null)
				this.writeFolder = FileUtils.removeExtention(this.expFile) + "/";

			init();
			this.filePathsExist();
			// normalize and center the data
			normalize();

			//make the covariance/correlation matrix
			String type = getCorrelationType();
			String covMatFN = this.writeFolder + "gene_" + type + ".txt";
			makeRelationShipMatrix(	covMatFN,
									type);

			// 1. calculate the correlation/covariance matrix and
			// 2. do the PCA,
			// 3. calculate cronbach alphas
			pca(type,
				covMatFN);
			break;
		}
		case 1://corrects PCs in a matrix from which the eigenvectors do not originate. Ensures same normalization (in rotate()) etc. is used.
		{
			// create output folder for PC corrected files
			checkAndMakeFolder(this.writeFolderCorrected);

			// put the matrix in the same space and calculate the PC scores for the PCs defined based on the public data
			MyMatrix[] rotationMatrixes = rotate();
			MyMatrix sampleStruct = rotationMatrixes[2];
			MyMatrix zScoreMatrix = null;
			//calculates Z-scores from the PCscores.

			String scoreFile = this.writeFolderCorrected + "SAMPLE.PC.scores.txt";

			MyMatrix scores = new MyMatrix(scoreFile);

			MyMatrix eigenVectors = rotationMatrixes[3];

			JuhaPCA.PCA.log("16. Adjusting for PCs");
			// for each of these numbers a file is created with that many PC's
			// corrected second number is replaced by this.PCs
			int[] PCAadjustments = new int[] {0,1,2,3,4,5,6,7,8,9,10, 50, 100, 200, 400,800};// ,132,300,500,1000};//,5000,eigenVectors.rows()};(2nd one gets replace by user input

			// if chr21 file is defined this will be used to determine the ranking of genes in this file vs those that are not this is determined after each PC correction iteratively
			adjustForPCs(	sampleStruct,
							PCAadjustments,
							eigenVectors,
							scores,
							this.writeFolderCorrected,
							this.writeFolder,
							zScoreMatrix,
							this.zScoresCutoff,
							this.log2,
							this.correctResultsForSTdevs,
							this.PCs,
							this.avgStdevFolder);

			System.out.println("Done, Results saved in: " + this.writeFolderCorrected);
			break;
		}
		case 2://correction gene for gene assuming PCA matrix. Works on matrixes with more then a billion entries as well. (unlike case 1:)
		{
			String centeredFN = writeFolder + "MATRIX_Centered.txt.gz";

			// create output folder for PC corrected files
			checkAndMakeFolder(this.writeFolderCorrected);

			File eigenVectorFn = null;
			if (this.isPcaOverGenes)
			{
				eigenVectorFn = new File(this.writeFolder + "GENE.eigenvectors.txt.gz");
				if (!eigenVectorFn.exists())
					eigenVectorFn = new File(this.writeFolder + "GENE.eigenvectors.txt");
			}
			else
			{
				eigenVectorFn = new File(this.writeFolder + "SAMPLE.eigenvectors.txt.gz");
				if (!eigenVectorFn.exists())
					eigenVectorFn = new File(this.writeFolder + "SAMPLE.eigenvectors.txt");
			}

			//read in matrix (genes are on the rows, PCs are on the columns in the eigenFile)
			ArrayList<Integer> pcList = parsePCs(PCs);
			p("userList=" + pcList);
			
			int largestPC = getLargest(pcList);
			p("largest PC corrected=" + largestPC);

			MatrixStruct eigenVectors = new MatrixStruct(	eigenVectorFn.getAbsolutePath(),
															largestPC + 1,
															-1);

			String scoreFn = this.writeFolderCorrected + "SAMPLE.PC.scores.txt";
			String preppedMatrixFN = centeredFN;

			preppedMatrixFN = FileUtils.removeExtention(preppedMatrixFN) + "_transposed.txt.gz";
			p("Transposing matrix");
			MyMatrix transposeMatrix = new MyMatrix(centeredFN);
			transposeMatrix.transpose();
			transposeMatrix.write(preppedMatrixFN);

			BufferedReader sampleReader = FileUtils.createReader(preppedMatrixFN);
			String header = sampleReader.readLine();//need to read the first line to put the filepointer past the headers
			String line = null;

			//create matrix that will hold a single sample
			MatrixStruct sampleStruct = new MatrixStruct(	preppedMatrixFN,
															1,
															-1);

			int gene = 0;

			p("Adjusting PCs:\t" + pcList);
			while ((line = sampleReader.readLine()) != null)
			{
				if (gene % 1000 == 0)
					p(gene + " genes finished");
				String[] rowCells = line.split("\t");
				String rowName = rowCells[0];
				double[] rowValues = new double[rowCells.length - 1];
				//parse row
				for (int c = 0; c < rowCells.length - 1; c++)
				{
					try
					{
						rowValues[c] = Double.parseDouble(rowCells[c + 1]);
					} catch (Exception e)
					{
						p("Unparsable value = " + rowValues[c]);
						e.printStackTrace();
					}
				}
				sampleStruct.setRow(0,
									rowName,
									rowValues);

				//calculate PC scores for this row
				MatrixStruct scores = PCA.scores1Row(	sampleStruct,
														eigenVectors,
														scoreFn,
														gene);

				for (int pc : pcList)// correct this sample for the selected PCs
				{
					if (pc > eigenVectors.rows())
						break;
					// remove the signal of this single PC from all the genes
					removeSignalAllEigenTransposed(	sampleStruct,
													scores,
													pc,
													gene,
													eigenVectors);
				}

				String writeFN = writeFolderCorrected + "PC_1-" + largestPC + ".txt.gz";
				if (this.getPCs() != null)
					writeFN = writeFolderCorrected + "PC_" + this.getPCs() + ".txt.gz";

				if (gene == 0)//create new file
				{
					p("Writing corrected data to:\t" + writeFN);
					sampleStruct.write(writeFN);
				}
				else//append
					sampleStruct.write(	writeFN,
										true);
				gene++;
			}

			//				//calculate andwrite cronbachs alphas
			//				p("Writing cronbach alphas at:\t" + cronbachFN);			
			//				String cronbachFN = scoreFn.toString().replace(".gz", "").replace(".txt", "") + ".cronbachsAlpha.txt";
			//				FileWriter fw = new FileWriter(cronbachFN);

			p("Done, Results saved in: " + this.writeFolderCorrected);
			break;
		}
		}
	}

	private void makeRelationShipMatrix(String covMatFN,
										String type)
	{
		if (!isPcaOverGenes)
			covMatFN = this.writeFolder + "sample_" + type + ".txt";

		//create the correlation matrix
		CorrelationLarge correlationLarge = new CorrelationLarge();
		correlationLarge.setExpressionFN(this.tempName);
		correlationLarge.setThreads(8);
		correlationLarge.setWriteFn(covMatFN);
		correlationLarge.setCorrelation(this.correlation);
		correlationLarge.run();
	}

	private String getCorrelationType()
	{
		String type = "covariance";
		if (this.correlation)
			type = "correlation";
		return type;
	}

	private MyMatrix[] rotate() throws IOException
	{
		/** 6. Calculate PCscores for single sample **/
		JuhaPCA.PCA.log(" 1. Loading sample matrix");
		MyMatrix singleSample = new MyMatrix(sampleFile);//expressionMatrix.getRow(0);
		singleSample.putGenesOnRows();

		//keep only the genes/rows that were used in the public samples as well
		String geneAveragesFN = writeFolder + "SAMPLE_Norm_GeneAverages.txt";
		MyMatrix geneAverages = new MyMatrix(geneAveragesFN);
		geneAverages.keepRows(singleSample);

		//rotate the sample to the same sample space
		singleSample = center(	singleSample,
								geneAverages,
								this);

		//calculate the PC scores
		//Get the eigenvector file based on the public data and put it in the right orientation.
		String saveNameSingleSampleScore = writeFolderCorrected + "SAMPLE.PC.txt";
		JuhaPCA.PCA.log(" 11. Loading gene eigen vector file: ");

		File eigenVectorFn = null;
		if (isPcaOverGenes)
		{
			eigenVectorFn = new File(writeFolder + "GENE.eigenvectors.txt.gz");
			if (!eigenVectorFn.exists())
				eigenVectorFn = new File(writeFolder + "GENE.eigenvectors.txt");
		}
		else
		{
			eigenVectorFn = new File(writeFolder + "SAMPLE.eigenvectors.txt.gz");
			if (!eigenVectorFn.exists())
				eigenVectorFn = new File(writeFolder + "SAMPLE.eigenvectors.txt");
		}

		//read in matrix (genes are on the rows, PCs are on the columns in the eigenFile)
		MyMatrix eigenVectors = new MyMatrix(	eigenVectorFn.getAbsolutePath(),
												-1,
												1002);

		if (eigenVectors.getRowHeaders()[0].contentEquals("PC1"))
			eigenVectors.transpose();

		if (isPcaOverGenes)
			geneAverages.keepRows(eigenVectors);//also changes the rows of geneEigenvectors to have the same positions as geneAverages

		singleSample = putInCorrectorrientation(singleSample,
												eigenVectors);

		if (!eigenVectors.getRowHeaders()[0].contentEquals("PC1"))
			eigenVectors.transpose();

		JuhaPCA.PCA.log(" 12. Calculate the PC scores: ");
		MyMatrix[] scoreResults = PCA.scores(	eigenVectors,
												singleSample,
												saveNameSingleSampleScore,
												false,
												false);

		JuhaPCA.PCA.log("Files Written to: " + writeFolderCorrected);
		return scoreResults;
	}

	private MyMatrix putInCorrectorrientation(	MyMatrix singleSample,
												MyMatrix eigenVectors)
	{
		if (isPcaOverGenes)
			singleSample.putGenesOnRows();
		else if (!hasAnyOverlap(singleSample.getRowHash(),
								eigenVectors.rowNames))
			singleSample.transpose();

		if (!hasAnyOverlap(	singleSample.getRowHash(),
							eigenVectors.rowNames))
		{
			System.out.println("Error: No resemblant rownames between eigenvector and sample rownames");
			System.exit(2);
		}
		return singleSample;
	}

	private static boolean hasAnyOverlap(	Hashtable<String, Integer> inputHash,
											String[] names)
	{
		for (String name : names)
		{
			if (inputHash.containsKey(name))
				return true;
		}
		return false;
	}

	private MyMatrix center(MyMatrix singleSample,
							MyMatrix geneAverages,
							PcaPipelineLite pcApipelineLite) throws IOException
	{
		makeFolder(writeFolderCorrected);
		JuhaPCA.PCA.log(" 3. Removing rows that do not exist in the averages vector from public data");
		geneAverages.keepRows(singleSample);

		if (deSeqNorm)
		{
			JuhaPCA.PCA.log(" 6. rLog transformation");
			MyMatrix geoMeans = new MyMatrix(writeFolder + "geoMean.txt");
			String swapFN = writeFolderCorrected + "swapFile.txt";
			singleSample.write(swapFN);
			DeSeqNorm.rLog(	singleSample,
							writeFolderCorrected,
							swapFN,
							geoMeans,
							null);
			singleSample.write(writeFolderCorrected + "rLogTransformed_" + deSeqNorm + ".txt");
		}

		if (log2)
		{
			//			if(correctTotalReadCount <= 0 && rLog <= 0) // need to add 1 before log to avoid log(0)
			//				addLogVal = 0.5;
			JuhaPCA.PCA.log(" 6. Log transforming");
			singleSample.log2Transform(addLogVal);//Doing this after the quantile normalization now
			singleSample.write(writeFolderCorrected + "normalized_log2.txt");
		}

		singleSample.write(writeFolderCorrected + "keepRows.txt");

		if (centerGenes)
		{
			JuhaPCA.PCA.log(" 9. Adjusting for gene averages (centering to target PC space)");
			singleSample.adjustForAverageAllGenes(geneAverages);
		}

		MyMatrix sampleAvgs = singleSample.getAveragesPerCol();
		String sampleAveragesFileName = writeFolderCorrected + "sampleAverages.txt";
		sampleAvgs.write(sampleAveragesFileName);
		if (centerSamples)
		{
			JuhaPCA.PCA.log(" 8. Adjusting for row averages (centering to target PC space)");
			singleSample.adjustForAverageAllSamples(sampleAvgs);
		}

		String centeredFN = writeFolderCorrected + "centered.txt";
		JuhaPCA.PCA.log(" 10. Writing PC centered file to: " + centeredFN);
		singleSample.write(centeredFN);
		return singleSample;
	}

	private int getLargest(ArrayList<Integer> userList)
	{
		int largest = 0;
		for (int number : userList)
		{
			if (number > largest)
				largest = number;
		}
		return largest;
	}

	private void checkAndMakeFolder(String writeFolderCorrected2)
	{
		if (sampleFile == null || sampleFile.length() == 0)
			sampleFile = expFile;
		if (this.writeFolderCorrected == null || this.writeFolderCorrected.length() == 0)
		{
			File sample = new File(this.sampleFile);
			this.writeFolderCorrected = this.getFolderName(this.writeFolder) + FileUtils.removeExtention(sample.getName()) + "/";
		}
		FileUtils.makeDir(this.writeFolderCorrected);
		p("Output folder for correced files:" + this.writeFolderCorrected);
	}

	private void selectGenes(	MyMatrix expressionStruct,
								String genesToInclude) throws IOException
	{
		if (genesToInclude != null)
			new MyMatrix(genesToInclude).keepRows(expressionStruct); // keeps them in the original order
	}

	public void normalize() throws IOException, NotConvergedException, InterruptedException
	{
		p(" 1. Reading expression file:\n" + this.expFile);
		MyMatrix expressionStruct = new MyMatrix(this.expFile);
		System.out.println("1=" + expressionStruct.rows());
		p("rows = " + expressionStruct.rows() + "cols = " + expressionStruct.cols());
		// transposes matrix if genes/transcripts are not on rows
		expressionStruct.putGenesOnRows();

		//remove rows without variance (problem in deseq normalization)
		expressionStruct.removeNoVariance();
		expressionStruct.write(this.writeFolder + "noVarianceRowsRemoved.txt.gz");
		
		// keep only a subset of genes
		selectGenes(expressionStruct,
					this.genesToInclude);

		// remove genes without variance
		//expressionStruct.removeNoVariance(this.writeFolder + "noVarRemoved.txt.gz");

		System.out.println("3=" + expressionStruct.rows());
		if (this.deSeqNorm)// Does not log the values, just does the DEseq based correction
		{
			// String writeGeoFN = this.writeFolder+ "geoMean.txt";
			DeSeqNorm.rLog(	this.writeFolder,
							expressionStruct,
							true,
							null);
		}

		if (this.log2)// log2 the data and add <this.addLogVal> before
		{
			LogTransform.log2(	this.writeFolder,
								expressionStruct,
								true,
								this.addLogVal);
		}

		JuhaPCA.PCA.log("12 Calculating STdevs");
		System.gc();
		MyMatrix stDevs = expressionStruct.stDevRows();
		stDevs.write(this.writeFolder + "gene_STDevs.txt");

		// expressionStruct.removeNoVariance(this.writeFolder+"noVarRemoved2.txt.gz");

		JuhaPCA.PCA.log("17. Calculating row averages");
		MyMatrix rowAverages = expressionStruct.getAveragesPerRow();
		rowAverages.write(this.writeFolder + "SAMPLE_Norm_GeneAverages.txt");
		if (this.centerGenes)
		{
			JuhaPCA.PCA.log("18. Centegering genes: Adjusting for gene averages");
			expressionStruct.adjustForAverageAllGenes(rowAverages);
		}

		JuhaPCA.PCA.log("15. Calculating column averages");
		MyMatrix colAverages = expressionStruct.getAveragesPerCol();
		String columnavgsfile = this.writeFolder + "SAMPLE_Norm_sampleAverages.txt";
		colAverages.write(columnavgsfile);
		if (this.centerSamples)
		{
			JuhaPCA.PCA.log("16. Centering samples: Adjusting for sample averages");
			expressionStruct.adjustForAverageAllSamples(colAverages);
			expressionStruct.write(this.writeFolder + "SAMPLE_adjustedForSampleAverages.txt.gz");
		}

		String expNormLogCentFile = this.writeFolder + "MATRIX_Centered.txt.gz";
		JuhaPCA.PCA.log("19. Writing centered file in: " + expNormLogCentFile);
		expressionStruct.write(expNormLogCentFile);

		if (!isPcaOverGenes)
			expressionStruct.putGenesOnCols();

		this.tempName = this.writeFolder + "pre_Correlation_Or_Covariance.txt";
		expressionStruct.write(this.tempName);
		expressionStruct = null;
		System.gc();
		System.gc();
		System.gc();
	}

	public void init() throws IOException, ParserConfigurationException, TransformerException
	{
		makeFolder(this.writeFolder);
		// create an XML file: //need to fix this so I can use it for calling Tessa's/Juha's prediction program that calculates AUCs per Reactome/GO term
		//		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		//		DocumentBuilder db = dbf.newDocumentBuilder();
		// create json File
		this.writeConfig();
	}

	void makeFolder(String writeFolder)
	{
		File folder = new File(writeFolder);
		if (!folder.exists())
		{
			folder.mkdir();
		}

	}

	void pca(	String type,
				String covMatFN) throws IOException, InterruptedException
	{
		System.gc();
		System.gc();
		System.gc();
		System.gc();

		Runtime runtime = Runtime.getRuntime();
		NumberFormat format = NumberFormat.getInstance();
		long maxMemory = runtime.maxMemory();
		long allocatedMemory = runtime.totalMemory();
		long freeMemory = runtime.freeMemory();

		System.out.println("free memory: " + format.format(freeMemory / 1024 / 1024 / 1024) + "<br/>");
		System.out.println("allocated memory: " + format.format(allocatedMemory / 1024 / 1024 / 1024) + "<br/>");
		System.out.println("max memory: " + format.format(maxMemory / 1024 / 1024 / 1024) + "<br/>");
		System.out.println("total free memory: " + format.format((freeMemory + (maxMemory - allocatedMemory) / 1024 / 1024 / 1024)));

		PCA.log("Matrix decomposition");
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(this.writeFolder + "export.sh")));

		writer.write("export DYLD_LIBRARY_PATH=\"/opt/intel/mkl/lib/\":\"/opt/intel/compilers_and_libraries_2016.1.111/mac/compiler/lib\":$DYLD_LIBRARY_PATH\n");
		writer.write("cd " + this.writeFolder + "\n");
		writer.write("/Volumes/Promise_RAID/juha/PCA++/pca evd " + covMatFN + "\n");
		writer.write("mv eigenvectors.txt " + this.writeFolder + "GENE.eigenvectors.txt" + "\n");
		writer.write("/Volumes/Promise_RAID/juha/PCA++/pca pc-scores " + type + " " + this.tempName + " GENE.eigenvectors.txt" + "\n");
		writer.write("mv eigenvalues.txt " + this.writeFolder + "GENE.eigenvalues.txt" + "\n");

		// get the number of PCs with Cronbach > 0.7
		String geneEigenVectorFN = "GENE.eigenvectors0.7.txt";
		writer.write("n=0\n" + "while read line; do \n" + "if [[ $line == \"Cronbach\"* ]]; then continue; fi\n" + "compare=$(echo $line'<'0.7 | bc)\n" + "if [[ compare -eq 1 ]]; then break; fi\n" + "((n=$n+1))\n" + "done < cronbach.txt\n" + "echo $n\n" + "cat < GENE.eigenvectors.txt | cut -f1-$n > " + geneEigenVectorFN + "\n");
		writer.write("gzip -f GENE.eigenvectors.txt\n");

		writer.close();
		String command = "sh " + this.writeFolder + "export.sh";
		runShell(command);

		if (new File(this.sampleFile).exists())
		{
			this.runMode = 1;
			this.run();
		}

		System.out.println("All done!");
	}

	private void runShell(String command)
	{
		System.out.println("Shellcommand = " + command);
		ExecCommand exec = new ExecCommand(command);
		System.out.println("execute output: \n" + exec.getOutput());
		System.out.println("execute error: \n" + exec.getError());
	}

	public String getWritePath(String name)
	{
		if (!name.contains("\\") && !name.contains("/"))
			return getFolderName(this.writeFolder) + name;
		return name;
	}

	public String getFolderName(String fn)
	{
		if (!fn.endsWith("/") && !fn.endsWith("\\"))
			fn = fn + "/";
		return fn;
	}

	public void readVars()
	{
		readVars(this.jsonFN);
	}

	private void adjustForPCs(	MyMatrix inputMatrix,
								int[] PCAadjustments,
								MyMatrix eigenVectors,
								MyMatrix scores,
								String writeFolder,
								String vectorFolder,
								MyMatrix zScores,
								double zScoresCutoff,
								boolean log2,
								boolean correctResultsForSTdevs,
								String PCs,
								String avgStdevFolder) throws IOException
	{
		if (isPcaOverGenes == false)
			inputMatrix.putGenesOnCols();
		// adjustOnZscores();
		MyMatrix tTestResults = null;
		MyMatrix difference = null;

		// get PCs to correct for from the user input if any
		ArrayList<Integer> userList = parsePCs(PCs);
		System.out.println("userList=" + userList);
		// JuhaPCA.PCA.log("Calculating variance explained by each PC"); first argument determins the number of principal components to
		// calculate the variance explained for. In this case it is the maximum number of principal components that is corrected for.
		// varianceExplained(PCAadjustments[PCAadjustments.length-1], writeFolder, scores, inputMatrix);//this function only works if the
		// gene standard deviation is set to 1 (and the average of course) JuhaPCA.PCA.log("Calculating variance explained done");
		for (int pcs = -1; pcs < PCAadjustments.length; pcs++)// for all the different numbers of PCs to correct for (e.g. int[] PCAadjustments = new int[]{0,100,2,25,300,500,1000};};)
		{
			MyMatrix sampleStruct = inputMatrix.copy();// make a copy so you can use the original for when you want to correct for a different number of PCs
			String writePCName = null;

			ArrayList<Integer> PCsToAdjust = new ArrayList<Integer>();
			if (pcs == -1)// if the user has a specific list he wants to correct
			{
				if (userList == null)
					continue;
				PCsToAdjust = userList;
				writePCName = PCs;
			}
			else
			{
				int adjustPCs = PCAadjustments[pcs];
				System.out.println("Adjusting:" + adjustPCs + " PCs");
				if (adjustPCs > eigenVectors.rows())// in case attempting to correct for more PCs than in the eigenvector file, just correct for all the eigenvectors in the eigenvector files.
					adjustPCs = eigenVectors.rows() + 1;
				writePCName = Integer.toString(adjustPCs);
				for (int p = 1; p < adjustPCs; p++)
				{
					PCsToAdjust.add(p);
				}
			}

			String writeFileName = writeFolder + "PC_1-" + writePCName + ".txt.gz";
			if (pcs != PCAadjustments.length - 1)// if it is not the last batch of PC corrections
			{
				correctPCs(	sampleStruct,
							scores,
							eigenVectors,
							PCsToAdjust);
			}

			System.out.println("writeFileName =" + writeFileName);
			sampleStruct.putGenesOnRows();
			sampleStruct.write(writeFileName);
			createExtraFiles(sampleStruct,
								correctResultsForSTdevs,
								vectorFolder,
								writeFileName);

		}
	}

	private void correctPCs(MyMatrix sampleStruct,
							MyMatrix scores,
							MyMatrix eigenVectors,
							ArrayList<Integer> PCsToAdjust)
	{
		p("PCsToAdjust = " + PCsToAdjust);
		for (int pc : PCsToAdjust)// correct this sample for the selected PCs
		{
			if (pc > eigenVectors.rows())
				break;
			for (int s = 0; s < sampleStruct.cols(); s++)
			{
				// remove the signal of this single PC from all the genes
				removeSignalAllgenes(	sampleStruct,
										scores,
										pc,
										s,
										eigenVectors,
										false);
			}
		}
	}

	private void removeSignalAllgenes(	MyMatrix sampleStruct,
										MyMatrix scores,
										int pc,
										int s,
										MyMatrix eigenVectors,
										boolean add)
	{
		for (int gene = 0; gene < sampleStruct.rows(); gene++)
		{
			double signal = scores.matrix.get(	pc - 1,
												s)
					* eigenVectors.matrix.get(	pc - 1,
												gene);
			sampleStruct.matrix.add(gene,
									s,
									-signal);
		}
	}

	private void removeSignalAll(	MatrixStruct sampleStruct,
									MatrixStruct scores,
									int pc,
									int gene,
									MatrixStruct eigenVectors)
	{
		// correct all the genes for this PC
		for (int sample = 0; sample < sampleStruct.cols(); sample++)
		{
			double signal = scores.matrix.get(	0,
												pc - 1)
					* eigenVectors.matrix.get(	sample,
												pc - 1);
			sampleStruct.matrix.add(0,
									sample,
									-signal);
		}

	}

	private void removeSignalAllEigenTransposed(MatrixStruct sampleStruct,
												MatrixStruct scores,
												int pc,
												int gene,
												MatrixStruct eigenVectors)
	{
		// correct all the genes for this PC
		for (int sample = 0; sample < sampleStruct.cols(); sample++)
		{
			double signal = scores.matrix.get(	0,
												pc - 1)
					* eigenVectors.matrix.get(	pc - 1,
												sample);
			sampleStruct.matrix.add(0,
									sample,
									-signal);
		}

	}

	private void createExtraFiles(	MyMatrix sampleStruct,
										boolean correctResultsForSTdevs,
										String vectorFolder,
										String writeFileName) throws IOException
	{
		// does not change the sampleStruct matrix
		devideBySTdev(	sampleStruct,
						correctResultsForSTdevs,
						vectorFolder,
						writeFileName);

		// adds averages and stdevs to matrix
//		String zscoreWriteFN = new File(writeFileName).getName().replace(	".txt",
//																			"_zScores.txt");
//		return Zscore.zScores(	this.writeFolderCorrected,
//								zscoreWriteFN,
//								sampleStruct,
//								this.avgStdevFolder,
//								zscoreWriteFN.replace(	"_zScores.txt",
//														"stats.txt"),
//								true);
		// int smoothNgenes = 10;
		// pca.PCA.log("18. Smooth "+smoothNgenes+" genes");
		// smoothSignal(sampleStruct.copy(),smoothNgenes, writeFileName);

		// smoothNgenes = 100;
		// JuhaPCA.PCA.log("19. Smooth "+smoothNgenes+" genes");
		// smoothSignal(sampleStruct.copy(),smoothNgenes,writeFileName);

		// double topPercent = 0.5;
		// keepTopPercentage(sampleStruct,
		// vectorFolder+"SAMPLE_Norm_columnAverages.txt",topPercent,
		// writeFileName, false);
		//
		// smoothNgenes = 10;
		// pca.PCA.log("21. Smooth "+smoothNgenes+" genes");
		// smoothSignal(sampleStruct.copy(),smoothNgenes,
		// writeFileName.replace(".txt", "Top"+topPercent*100+"%_.txt"));
		//
		// smoothNgenes = 100;
		// pca.PCA.log("22. Smooth "+smoothNgenes+" genes");
		// smoothSignal(sampleStruct.copy(),smoothNgenes,writeFileName.replace(".txt",
		// "Top"+topPercent*100+"%_.txt"));

	}

	void keepTopPercentage(	MyMatrix sampleStruct,
							String averagesFN,
							double topPercent,
							String writeFileName,
							boolean lowest,
							boolean writeAll) throws IOException
	{
		JuhaPCA.PCA.log("20. Highest " + topPercent * 100 + "% only");
		MyMatrix averages = new MyMatrix(averagesFN);
		averages.sortCol(0);
		averages.write(averagesFN.replace(	".txt",
											"SAMPLE_Norm_columnAverages_SORTED_ALL.txt"));
		MyMatrix part = keepPart(	averages,
									topPercent,
									lowest);
		averages.write(averagesFN.replace(	".txt",
											"SAMPLE_Norm_columnAverages_SORTED.txt"));
		sampleStruct.keepRows(averages);
		part.write(averagesFN.replace(	".txt",
										"SAMPLE_Norm_columnAverages_" + topPercent + "highestExpressed.txt"));
		sampleStruct.keepRows(part);
		if (writeAll)
			sampleStruct.write(writeFileName.replace(	".txt",
														".Top" + topPercent * 100 + "%only.txt"));
	}

	private void devideBySTdev(	MyMatrix sampleStruct,
								boolean correctResultsForSTdevs,
								String vectorFolder,
								String writeFileName) throws IOException
	{
		MyMatrix copy = sampleStruct.copy();

		if (correctResultsForSTdevs)
		{
			JuhaPCA.PCA.log("16. Divide by standard deviation");
			MyMatrix STdevs = new MyMatrix(vectorFolder + "gene_STDevs.txt");
			STdevs.keepRows(copy);
			// if the standard deviation is smaller than 1, set it to 1 to avoid inflated values for genes that have a very small stdev
			for (int s = 0; s < STdevs.rows(); s++)
				if (STdevs.matrix.get(	s,
										0) < 1)
					STdevs.matrix.set(	s,
										0,
										1);
			copy.divideBy(	STdevs,
							true);
			copy.write(writeFileName.replace(	".txt",
												"DevidedBySTdevs.txt"));
		}

	}

	private MyMatrix keepPart(	MyMatrix averages,
								double topPercent,
								boolean lowest) // returns the last part of the matrix
	{
		int topX = (int) (topPercent * averages.rows());
		MyMatrix part = new MyMatrix(	topX,
										averages.cols());
		part.setColHeaders(averages.getColHeaders());
		if (!lowest)
			for (int r = 0; r < topX; r++)
			{
				part.setRow(r,
							averages.getRowHeaders()[r],
							averages.getRowValues(r));
			}
		else
			for (int r = averages.rows() - topX, out = 0; r < averages.rows(); r++, out++)
			{
				part.setRow(out,
							averages.getRowHeaders()[r],
							averages.getRowValues(r));
			}
		return part;
	}

	private void smoothSignal(	MyMatrix sampleStruct,
								int n,
								String writeFileName) throws IOException
	{
		if (n > sampleStruct.rows())
			return;
		if (n % 2 == 0)// want an uneven number so you can take half before half after + the one itself
		{
			n++;
			// System.out.println("Taking average of " + n + " numbers");
		}
		for (int c = 0; c < sampleStruct.cols(); c++)
		{
			double runningSum = 0;
			double runningAvg = 0;
			double[] avgs = new double[sampleStruct.rows()];
			int count = 0;

			for (int r = 0; r < sampleStruct.rows(); r++)
			{
				if (count < n)// for the first n numbers
				{
					runningSum += sampleStruct.matrix.get(	r,
															c);
					count++;
				}
				else
				{
					runningSum -= sampleStruct.matrix.get(	r - n,
															c);
					runningSum += sampleStruct.matrix.get(	r,
															c);
				}
				runningAvg = runningSum / count;
				if (r - n / 2 >= 0)
					avgs[r - n / 2] = runningAvg;
			}
			// last n/2 rows
			for (int r = sampleStruct.rows(); r < sampleStruct.rows() + n / 2; r++)
			{
				runningSum -= sampleStruct.matrix.get(	r - n,
														c);
				count--;
				runningAvg = runningSum / count;
				avgs[r - n / 2] = runningAvg;
			}
			sampleStruct.setCol(avgs,
								c);
		}
		sampleStruct.write(writeFileName.replace(	".txt",
													".Smoothed" + n + "Genes.txt"));
	}

	public ArrayList<Integer> parsePCs(String PCsToAdjust)
	{// function takes format like: "1,4,5-10,3-10"
		ArrayList<Integer> PCs = new ArrayList<Integer>();
		if (PCsToAdjust == null || PCsToAdjust.contains("null"))
		{
			for (int i = 1; i <= 100; i++)
				PCs.add(i);
			return PCs;
		}
		
		if (PCsToAdjust.contains("-"))
		{
			String[] eles = PCsToAdjust.split("-");
			for (int e = 0; e < eles.length - 1; e++)
			{
				String[] ele = eles[e].split(",");
				int last = ele.length - 1;
				int start = Integer.parseInt(ele[last]);
				String[] ele2 = eles[e + 1].split(",");
				int end = Integer.parseInt(ele2[0]);
				for (int n = start; n <= end; n++)
					PCs.add(n);
			}
		}
		if (PCsToAdjust.contains(","))
		{
			String[] eles = PCsToAdjust.split(",");
			for (int e = 0; e < eles.length; e++)
			{
				if (!eles[e].contains("-"))
					PCs.add(Integer.parseInt(eles[e]));
			}
		}

		// System.out.println("size = " + PCs.size());
		// for(int PC : PCs)
		// {
		// System.out.println("pc = " + PC);
		// }

		return PCs;
	}

	// Prediction script
	public boolean filePathsExist()
	{
		File[] files = new File[1];

		files[0] = new File(this.expFile);

		for (File file : files)
		{
			if (file == null)
			{
				System.out.println("NOT ALL NECESSARY FILES HAVE BEEN SUPPLIED!\n");
				return false;
			}
			if (file.exists())
				continue;
			System.out.println("THIS FILE/FOLDER DOES NOT EXIST!\n" + file.getAbsolutePath());
			return false;
		}
		return true;
	}

	public int getRunMode()
	{
		return runMode;
	}

	public void setRunMode(int runMode)
	{
		this.runMode = runMode;
	}

	public String getExpFile()
	{
		return expFile;
	}

	public void setExpFile(String expFile)
	{
		this.expFile = expFile;
	}

	public String getWriteFolder()
	{
		return writeFolder;
	}

	public void setWriteFolder(String writeFolder)
	{
		this.writeFolder = writeFolder;
	}

	public String getGenesToInclude()
	{
		return genesToInclude;
	}

	public void setGenesToInclude(String genesToInclude)
	{
		this.genesToInclude = genesToInclude;
	}

	public boolean isPcaOverGenes()
	{
		return isPcaOverGenes;
	}

	public void setPcaOverGenes(boolean isPcaOverGenes)
	{
		this.isPcaOverGenes = isPcaOverGenes;
	}

	public boolean isCorrelation()
	{
		return correlation;
	}

	public void setCorrelation(boolean correlation)
	{
		this.correlation = correlation;
	}

	public boolean isCenterGenes()
	{
		return centerGenes;
	}

	public void setCenterGenes(boolean centerGenes)
	{
		this.centerGenes = centerGenes;
	}

	public boolean isCenterSamples()
	{
		return centerSamples;
	}

	public void setCenterSamples(boolean centerSamples)
	{
		this.centerSamples = centerSamples;
	}

	public boolean isLog2()
	{
		return log2;
	}

	public void setLog2(boolean log2)
	{
		this.log2 = log2;
	}

	public double getAddLogVal()
	{
		return addLogVal;
	}

	public void setAddLogVal(double addLogVal)
	{
		this.addLogVal = addLogVal;
	}

	public boolean isDeSeqNorm()
	{
		return deSeqNorm;
	}

	public void setDeSeqNorm(boolean deSeqNorm)
	{
		this.deSeqNorm = deSeqNorm;
	}

	public String getSampleFile()
	{
		return sampleFile;
	}

	public void setSampleFile(String sampleFile)
	{
		this.sampleFile = sampleFile;
	}

	public String getPCs()
	{
		return PCs;
	}

	public void setPCs(String pCs)
	{
		PCs = pCs;
	}

	public String getWriteFolderCorrected()
	{
		return writeFolderCorrected;
	}

	public void setWriteFolderCorrected(String writeFolderCorrected)
	{
		this.writeFolderCorrected = writeFolderCorrected;
	}

	public String getAvgStdevFolder()
	{
		return avgStdevFolder;
	}

	public void setAvgStdevFolder(String avgStdevFolder)
	{
		this.avgStdevFolder = avgStdevFolder;
	}

	public double getzScoresCutoff()
	{
		return zScoresCutoff;
	}

	public void setzScoresCutoff(double zScoresCutoff)
	{
		this.zScoresCutoff = zScoresCutoff;
	}

	public boolean isCorrectResultsForSTdevs()
	{
		return correctResultsForSTdevs;
	}

	public void setCorrectResultsForSTdevs(boolean correctResultsForSTdevs)
	{
		this.correctResultsForSTdevs = correctResultsForSTdevs;
	}

	public String getTempName()
	{
		return tempName;
	}

	public void setTempName(String tempName)
	{
		this.tempName = tempName;
	}

	public static long getSerialversionuid()
	{
		return serialVersionUID;
	}

}

// private static Matrix inDirectPCA(Matrix expressionStruct,
// String writeFolder, boolean correlation, double[] cutoffs) throws
// IOException, NotConvergedException {
// if(expressionStruct.getRowHeaders()[0].contains("ENSG0") &&
// expressionStruct.getRowHeaders()[0].contains("ENST0"))
// {
// System.out.println("Genes should be on Cols and are not, transposing (assumes
// geneIDs start with ENSG0 or ENST0)");
// expressionStruct.transpose();
// }
// String type = "covariance";
//// expressionStruct.removeNoVariance(writeFolder+"noVarRemoved.txt");
//
// if(correlation)
// {
// type = "correlation";
// }
// ConcurrentCovariation calculator = new ConcurrentCovariation(20);
// double[][] inMat = expressionStruct.getMatrix();
// double[][] covMatrix = calculator.pairwiseCovariation(inMat,false, null,
// expressionStruct.getRowHeaders(),correlation, cutoffs);
//
// Matrix covMat = new Matrix(expressionStruct.getRowHeaders(),
// expressionStruct.getRowHeaders(), covMatrix);
//
//
// pca.PCA.log("21. Writing covariance matrix over the samples");
// String covMatFN = writeFolder+"SAMPLE_"+type+".txt.gz";
// covMat.write(covMatFN);
//
//// to test something
// Matrix covMatOtherFormat = new Matrix(expressionStruct.getRowHeaders(),
// expressionStruct.getRowHeaders(), covMatrix);
// covMatOtherFormat.getAverageCols(true)
// .write(covMatFN.replace(".txt", "_Absolute_averages.txt"));
// covMatOtherFormat.write(covMatFN);
//
// inMat = null; covMatrix = null; System.gc();System.gc();
//
// pca.PCA.log("22. calculating eigenvalues over the samples");
// Matrix[] evds = pca.PCA.evd(covMat, Paths.get(writeFolder+"SAMPLE"));
// Matrix eigenVectors = evds[0];
// Matrix PCeigenvalues = evds[1];
//
// pca.PCA.log("23. calculating PCscores over the genes");
// String saveNamePCscoresGene = writeFolder + "GENE_PC.txt";
// Matrix[] PCscoresGenesAndAverages =
// PCA.scores(eigenVectors,expressionStruct, saveNamePCscoresGene);
// Matrix PCscoresGenes = PCscoresGenesAndAverages[0];
// PCscoresGenesAndAverages=null;
// System.gc();
//
// pca.PCA.log("24. Transform PCscores to eigenvectors of Genes");
// String saveNameEigenVectorsOverGenes = writeFolder +
// "GENE.eigenvectors.txt.gz";
// Matrix geneEigenVectors = PCA.transform(PCscoresGenes, PCeigenvalues,
// saveNameEigenVectorsOverGenes);
// return geneEigenVectors;
// }
