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
import Tools.ExecCommand;
import Tools.Script;
import eqtlmappingpipeline.normalization.Normalizer;
import no.uib.cipr.matrix.NotConvergedException;

public class PCApipeline extends Script<PCApipeline>
{
	/*
	 * 1. load expression matrix //2. removeDuplicates (any samples with correlation above <removeDuplicates> threshold) 3. put genes on rows if they are not already 4. Sort genes by chromosome and remove any that are not on in the <chromLocationsFile> //5. add random values (add a random value of <randomAddValue> to any gene with an expression lower than <randomAddValue>) //6. remove low expression (remove genes that do not have an expression of at least <minExpression> in at least <minSamplesExpressed> samples. //7. remove a single gene with name <gene> (if you want to) //8. quantile normalize //8. Voom normalize 8. Rlog normalization (without logging data). This is referring to part of the normalization DEseq uses 9. 2log transform (adds <value> before doing transformation) //10. Take the <highestExpressed> % genes //12. keep genes with highest <TopVariance> % variance //13. correct Input for standard deviation if <correctInputForSTdevs> == true //14. calculate spearman correlation if <spearman> >=0. Uses <spearman> as minimum expression (if lower, expression is set to this number) //15. if <setLowestToAverage> == True, all the lowest values to the average (if allowing log(0), all 0 will become this average)) 16. Calculate column average and store them 17. Correct for column averages 18. Calculate rowaverages and store them //19. Correct for row averages //20. correctInputForSTdevsAfterCenter
	 * 
	 * 21. Put genes on rows if they are not already 22. Remove genes that have no variance 23. calculate correlation/covariance matrix 24. PCA over matrix 25. Calculate cronbach alphas (and therefore PCscores over samples) Different script(s): 26. Take all PCs>0.7 and recalculate gene eigenvectors 27. Put in database
	 */
	// static CreateGeneEigenvectorFileVariables Var = new CreateGeneEigenvectorFileVariables();
	private static final long serialVersionUID = -8656788330691209748L;
	int runMode = 0; // 0 = PCA, 1 = datacorrection
	public String expFile = "/Volumes/Promise_RAID/GeneNetwork/Sipko/TestExpression2/TESTexpression.txt"; // expression file
	public String chromLocationsFile = null; // File that contains the chromosome locations
	public String writeFolder = "E:/Groningen/Data/Juha/Genes31995/31.07.pc1.illumina.genes.expressed.DEseqnorm/PCcorrection/"; // Folder where to write
	public String xmlFN = "config.xml"; // don't touch this
	public String GCgenes = null; // Not tested, don't use. If a file is supplied with the GC content per gene it will correct for this
	public String removeGene = null; // if there is a particular gene you want to remove from the matrix for any reason
	public String genesToInclude = null; // List of genes. This list is selected prior to any other steps. Other genes are discarded

	public boolean isPcaOverGenes = true; // if false it does the PCA over the samples rather than the genes
	public boolean correlation = false; // if false uses covariance
	public boolean setLowestToAverage = false; // sets all the lowest values in a sample to the average, effectively this means any gene that has an expression of does not contribute to the covariance or correlation
	public boolean centerSamples = false; // adjusts the expression of genes so that the average expression of a sample become 0 (centering over the samples)
	public boolean centerGenes = true; // adjusts the expression of genes so that the average expression of a sample become 0 (centering over the samples)
	public boolean writeAll = true; // write all intermediate files, is slower but helps finding understanding what happens in each step
	public boolean correctInputForSTdevs = false; // corrects the input for the standard deviation, can be done over genes or samples (look at the function itself it has a true/false argument)
	public boolean correctInputForSTdevsAfterCenter = false; // same as previous, but this time after centering the data
	public boolean log2 = true;
	public boolean skipQuantileNorm = true;
	public boolean STdevCutoff = false; // calculates the standard deviation of all genes and throws out the <highestExpressed> percent genes with the lowest standard deviation (instead of using average expression for cutoff)
	public boolean zScores = false;
	public boolean directPCA = true; // PCA over genes rather than samples.
	public boolean rLog = true; // If true uses DEseq normalization

	public int minSamplesExpressed = -1; // if -1, this will include all samples. Otherwise will exclude any gene that is expressed in less then this number of samples (and those that have no variance) any gene that has an expression lower then "minExpression" is considered as "not expressed"
	public int minExpression = -1; // if -1 does nothing. Otherwise sets the cutoff for minSamplesExpressed; see "minSamplesExpressed" public int ignoreLowestValues = -1;//Sets lowest values to the average so it does not contribute toward positive correlation.

	public double addLogVal = 0.5; // Value to add before taking the logarithm
	public double correctTotalReadCount = 0; // log((gene+0.5)/total*value)

	public double randomValue = 0; // if 0 does nothing. Otherwise adds a random value that is below this value, to values below this value.
	public double duplicateCutoff = 1; // if 1 does nothing. Samples with a correlation above this value are removed (only samples next to each other are removed, as this is usually the case for replicates, which we aim to remove like this, but not others)
	public double highestExpressed = 1; // 1 = all genes, 0.5 = 50% highest expressed genes only (removes 50% lowest expressed genes after quantile normalization (then re-normalizes)).
	public double spearman = -1; // if 0, does spearman, if >0, it sets all values below this value to 0.

	public String correlationScript = "/Volumes/Promise_RAID/GeneNetwork/Sipko/CorrelationLargeGenes.jar"; // script used to calculate correlation.

	//Prediction script
	public String geneNameFile = "/Volumes/Promise_RAID/GeneNetwork/ENSGToGeneNameHGNCBiotypeChromosomeStartStopStrandAndDescriptionV75.txt.filtered.txt";
	public String itemSetFile = "/Volumes/Promise_RAID/GeneNetwork/Sipko/Scripts/HPO/HPO.gmt";
	public String geneTermOutFile = "HPO_geneterm.txt";
	public String termDesCoutFile = "HPO_terms.txt";
	public String itemType = "Gene";
	public String limitToItemsInGenesetFile = "false"; // Data/OldSchoolExpression/GeneAnnotation/Genes.txt
	public String maxItems = "10000"; // Maximum number of items allowed in Tessa script
	public String minItems = "10"; // Not really sure how this is different from absMinItems
	public String absMinItems = "5"; // absolute minimum number of items in geneterm
	public String useTtest = "true"; // if false Wilcoxon will be used
	public String runRealAnalysis = "true";
	public String nrPermutations = "0";
	public String writeBinary = "false";
	public String logToFile = "false";
	public String label = "HPO-unscaled";
	public String javaExe = "/usr/bin/java -jar";
	public String RNAseqJar = "/Volumes/Promise_RAID/juha/PublicRNASeq/Prediction/RNASeq.jar";
	public String populateGenesetDBjs = "/Volumes/Promise_RAID/GeneNetwork/Sipko/05-2016/Tessa/populateGenesetDBTXT.js";
	public String geneDB = "genedb";
	public String hpoDB = "hpodb";
	public String termtype = "HPO";
	public String tessaPredictionFolder = "/Volumes/Promise_RAID/GeneNetwork/Sipko/05-2016/Tessa/Server_testing_DB/";
	public String tessaPredictionScript = "Ranking_multipleVersionsCombined.js";
	public String tessaMaxGenesPerTerm = "max15";
	public String tessaType = "regular";
	// if these variables are set the PC correction will be run after the eigenvectors are created (or if eigenvectors already exist only this part is run)
	public String sampleFile = "";
	public String chr21FN = "";
	public String PCs = "";// function takes format like: "1,4,5-10,3-10"
	public String writeFolderCorrected = "";
	public String avgStdevFolder = null; // this is the folder with the files for
	// the z-score calculations containing the averages and standard deviations to be used for this. allows z-scores to be calculated using the BBMRI averages and stdevs
	public double zScoresCutoff = Double.parseDouble("0");
	public boolean correctResultsForSTdevs = true;

	public int optimalPCremoval = -1;
	public String tempName;

	// System.out.println("This script calculates the eigenvectors over the
	// genes and uses the following input:\n"
	// + "fileName=<fileName> - Expression file (samples on rows, gene names on
	// columns)\n"
	// + "chrom=<fileName chromosome file> - File with genes names in 1st
	// column, chromosome in 2nd, position in 3rd (default=null)\n"
	// + "writeAll=<false/true> - writes all files from intermediary steps
	// (default=true)\n"
	// + "correctSTdevs=<false/true> - Corrects for STdevs prior to calculating
	// covariance matrix (default=false)\n"
	// + "log2=<false/true> - Log2 transformation after quantile normalization
	// (default=true)\n"
	// + "randomValue=<number> - adds a random value*<number> to numbers smaller
	// then <number> (default=1)\n"
	// + "writeFolder=<folderName> - folder to write the results in
	// (default=filename.replace(.txt,/)\n"
	// + "duplicate=<number> - correlation cutoff to consider samples as
	// duplicates (default=0.999)\n"
	// + "skipQN=<false/true> - whether to skip the Quantile normalizatin or not
	// (default=true)\n"
	// + "DESeqNorm=<false/true> - whether to use DESeq normalization
	// (default=true)\n"
	// + "highestExpressed=<number> - Percentage of genes to keep (genes most
	// highly expressed on average)"
	// + "genestoinclude=<filename> - File with in the first column the genes to
	// keep (assumes column headers)");

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
		this.writeConfig(this.jsonFN,this, false,true);
		switch (runMode)
		{
		case 0:
		{
			if (this.writeFolder == null)
				this.writeFolder = this.expFile.replace(".txt",
														"").replace(".gz",
																	"")
						+ "/";

			writeParameters();
			this.filePathsExist();
			// normalize and center the data
			normalize();
			// 1. calculate the correlation/covariance matrix and
			// 2. do the PCA,
			// 3. calculate cronbach alphas
			pca();
			break;
		}
		case 1:
		{
			// create output folder for PC corrected files
			if (this.writeFolderCorrected == null || this.writeFolderCorrected.length() == 0)
			{
				File sample = new File(this.sampleFile);
				this.writeFolderCorrected = this.getFolderName(this.writeFolder) + sample.getName().replace(".txt",
																											"").replace(".gz",
																														"")
						+ "/";
			}
			makeFolder(this.writeFolderCorrected);

			// put the matrix in the same space and calculate the PC scores for the PCs defined based on the public data
			MatrixStruct[] rotationMatrixes = RotateSample.rotate(this);
			MatrixStruct sampleStruct = rotationMatrixes[2];
			System.out.println("rows = " + sampleStruct.rows());

			MatrixStruct zScoreMatrix = null;

			// converts the PCscores to Z-scores. Bit slow and not really using it at the moment so commented it out.
			// pca.PCA.log("13. Calculating zScores");
			// String zScoreStats = this.writeFolder+"pcZscores_Stats.txt";
			// zScoreMatrix = rotationMatrixes[0].copy();
			// Zscore.changeToZscores(zScoreMatrix, zScoreStats);
			//
			// pca.PCA.log("14. Writing zScores");
			// zScoreMatrix.write(this.writeFolderCorrected+"pcZscoresSamples.txt");

			String scoreFile = this.writeFolderCorrected + "SAMPLE.PC.scores.txt";
			MatrixStruct scores = new MatrixStruct(scoreFile);

			MatrixStruct eigenVectors = rotationMatrixes[3];

			JuhaPCA.PCA.log("16. Adjusting for PCs");
			// for each of these numbers a file is created with that many PC's
			// corrected second number is replaced by this.PCs
			int[] PCAadjustments = new int[] { 0, 2, 100, 300};// ,132,300,500,1000};//,5000,eigenVectors.rows()};(2nd one gets replace by user input

			// if chr21 file is defined this will be used to determine the ranking of genes in this file vs those that are not this is determined after each PC correction iteratively
			MatrixStruct chr21 = null;
			if (this.chr21FN != null && new File(this.chr21FN).exists())
				chr21 = new MatrixStruct(this.chr21FN);

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
							chr21,
							this.optimalPCremoval,
							this.PCs,
							this.avgStdevFolder);

			System.out.println("Done, Results saved in: " + this.writeFolderCorrected);
		}

		}
	}

	private void selectGenes(	MatrixStruct expressionStruct,
								String genesToInclude) throws IOException
	{
		if (genesToInclude != null)
			new MatrixStruct(genesToInclude).keepRows(expressionStruct); // keeps them in the original order
	}

	public void normalize() throws IOException, NotConvergedException, InterruptedException
	{
		JuhaPCA.PCA.log(" 1. Reading expression file");
		MatrixStruct expressionStruct = new MatrixStruct(this.expFile);
		// transposes matrix if genes/transcripts are not on rows

		expressionStruct.putGenesOnRows();

		// keep only a subset of genes
		selectGenes(expressionStruct,
					this.genesToInclude);

		if (this.writeAll)
			expressionStruct.write(this.writeFolder + "startMatrix.txt.gz");

		// removes duplicates if this.duplicateCutoff <1
		RemoveDuplicates.removeDuplicates(	expressionStruct,
											this.duplicateCutoff,
											this.writeFolder,
											this.writeAll);

		// sorting on chromosome locations
		expressionStruct = SortChromosome.sort(	expressionStruct,
												this.chromLocationsFile);

		// adding random values, if this.randomValue =< 0, skips this
		expressionStruct.addRandomValues(	this.randomValue,
											this);

		System.out.println("Rows1 = " + expressionStruct.rows() + " cols = " + expressionStruct.cols());
		// remove genes without variance

		expressionStruct.removeNoVariance(this.writeFolder + "noVarRemoved.txt.gz");
		// remove genes that have an expression below this.minExpression
		if (this.minExpression > 0)// if this.minExpression is zero or smaller, skip this
			expressionStruct.removeLowExpression(	(this.writeFolder + "RemovedBelow" + this.minExpression + "inAtLeast" + this.minSamplesExpressed + "samples.txt.gz"),
													this.minSamplesExpressed,
													this.minExpression);
		System.out.println("Rows2 = " + expressionStruct.rows() + " cols = " + expressionStruct.cols());
		if (this.removeGene != null)// remove the gene defined by this.removeGene
			expressionStruct.removeRow(this.removeGene);

		if (!this.skipQuantileNorm)// quantile Normalize if false
			QuantileNormalize.quantileNormalize(expressionStruct,
												this.writeFolder,
												this.writeAll);

		if (this.correctTotalReadCount > 0)// correct for total number of reads
			CorrectReadcounts.correct(	this.writeFolder,
										this.correctTotalReadCount,
										expressionStruct,
										this.writeAll,
										0.5);

		if (this.rLog)// Does not log the values, just does the DEseq based correction
		{
			// String writeGeoFN = this.writeFolder+ "geoMean.txt";
			RLog.rLog(	this.writeFolder,
						expressionStruct,
						this.writeAll,
						null);
		}

		if (this.log2)// log2 the data and add <this.addLogVal> before
			LogTransform.log2(	this.writeFolder,
								expressionStruct,
								this.writeAll,
								this.addLogVal);

		if (this.highestExpressed > 0 && this.highestExpressed < 1)// keep the highest <this.highestExpressed> genes only
			HighestExpressed.highestExpressed(	expressionStruct,
												this.skipQuantileNorm,
												this.correctTotalReadCount,
												this.writeFolder,
												this.highestExpressed,
												this.STdevCutoff,
												this.writeAll);

		JuhaPCA.PCA.log("12 Calculating STdevs");
		System.gc();
		System.gc();
		MatrixStruct stDevs = expressionStruct.stDevRows();
		stDevs.write(this.writeFolder + "gene_STDevs.txt");
		System.out.println("Rows = " + expressionStruct.rows() + " cols = " + expressionStruct.cols());

		if (this.correctInputForSTdevs)
		{
			JuhaPCA.PCA.log("13 Divide all gene values by STdev for each gene ");
			expressionStruct.divideBy(	stDevs,
										false);// false corrects columns, true corrects rows

			JuhaPCA.PCA.log("14 Writing matrix divided by gene STdevs");
			expressionStruct.write(this.writeFolder + "_DividedBySGenesSTdev.txt.gz");
		}

		// expressionStruct.isExpressed(null, expressionStruct.cols()*0.05,10);//genes need to have an expression of at least 1 or more in at least 5% of all the samples.
		if (this.spearman >= 0)
		{
			Spearman.ranks(	this.writeFolder,
							expressionStruct,
							this.spearman);
			this.correlation = true;
		}

		if (this.setLowestToAverage)// this will cause the lowest values not to contribute to correlation or covariance
			LowestToAverage.lowestToAverage(expressionStruct);
		if (this.GCgenes != null)
		{
			MatrixStruct geneGC = new MatrixStruct(this.GCgenes);
			expressionStruct = GCcontent.calculateAndCorrect(	expressionStruct,
																geneGC,
																this.writeFolder + "gCperSampleWriteFN.txt",
																this.writeFolder + "GCcorrected.txt.gz");
		}

		// expressionStruct.removeNoVariance(this.writeFolder+"noVarRemoved2.txt.gz");

		JuhaPCA.PCA.log("15. Calculating column averages");
		MatrixStruct colAverages = expressionStruct.getAveragesPerCol();
		String columnavgsfile = this.writeFolder + "SAMPLE_Norm_sampleAverages.txt";
		colAverages.write(columnavgsfile);
		if (this.centerSamples)
		{
			JuhaPCA.PCA.log("16. Centering samples: Adjusting for sample averages");
			expressionStruct.adjustForAverageAllSamples(colAverages);
			expressionStruct.write(this.writeFolder + "SAMPLE_adjustedForSampleAverages.txt.gz");
		}
		JuhaPCA.PCA.log("17. Calculating row averages");
		MatrixStruct rowAverages = expressionStruct.getAveragesPerRow();
		rowAverages.write(this.writeFolder + "SAMPLE_Norm_GeneAverages.txt");
		if (this.centerGenes)
		{
			JuhaPCA.PCA.log("18. Centegering genes: Adjusting for gene averages");
			expressionStruct.adjustForAverageAllGenes(rowAverages);
		}

		String expNormLogCentFile = this.writeFolder + "MATRIX_Centered.txt.gz";
		JuhaPCA.PCA.log("19. Writing centered file in: " + expNormLogCentFile);
		expressionStruct.write(expNormLogCentFile);

		if (this.correctInputForSTdevsAfterCenter)
		{
			MatrixStruct stDevsRows = expressionStruct.stDevRows();
			stDevsRows.write(this.writeFolder + "_SampleStDevs.txt");
			JuhaPCA.PCA.log("13 Divide all gene values by STdev for each sample");
			expressionStruct.divideBy(	stDevsRows,
										true);// false corrects columns, true corrects rows

			JuhaPCA.PCA.log("14 Writing matrix divided by gene STdevs");
			expressionStruct.write(this.writeFolder + "_DividedBySGenesSTdev.txt.gz");
		}

		this.tempName = this.writeFolder + "pre_Correlation_Or_Covariance.txt";
		if(!isPcaOverGenes)
			expressionStruct.putGenesOnCols();
		expressionStruct.write(this.tempName);
		expressionStruct = null;
		System.gc();
		System.gc();
		System.gc();
	}

	public void writeParameters() throws IOException, ParserConfigurationException, TransformerException
	{
		makeFolder(this.writeFolder);
		DateFormat dateFormat = new SimpleDateFormat("dd/MM/yyyy HH:mm:ss");
		Date date = new Date();
		// create an XML file: //need to fix this so I can use it for calling Tessa's/Juha's prediction program that calculates AUCs per Reactome/GO term
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		DocumentBuilder db = dbf.newDocumentBuilder();
		// create json File
		//this.jsonFN = this.writeFolder + "config.json";	
		this.writeConfig();
		// if(this.jsonFN == null)
		// System.exit(1);
	}

	private void writeXML()
	{
//		String writeFN = this.writeFolder + getXML("xmlfile");
//		TransformerFactory tff = TransformerFactory.newInstance();
//		Transformer tf;
//		try
//		{
//			tf = tff.newTransformer();
//
//			DOMSource source = new DOMSource(xml);
//			StreamResult result = new StreamResult(new File(writeFN).getAbsolutePath());
//			tf.setOutputProperty(	OutputKeys.INDENT,
//									"yes");
//			tf.setOutputProperty(	"{http://xml.apache.org/xslt}indent-amount",
//									"2");
//			tf.transform(	source,
//							result);
//		} catch (TransformerConfigurationException e)
//		{
//			e.printStackTrace();
//		} catch (TransformerException e)
//		{
//			e.printStackTrace();
//		}
	}

	void makeFolder(String writeFolder)
	{
		File folder = new File(writeFolder);
		if (!folder.exists())
		{
			folder.mkdir();
		}

	}

	void pca() throws IOException, InterruptedException
	{
		System.gc();
		System.gc();
		System.gc();
		System.gc();

		String type = "covariance";
		if (this.correlation)
			type = "correlation";
		String covMatFN = this.writeFolder + "gene_" + type + ".txt";
		if(!isPcaOverGenes)
			covMatFN = this.writeFolder + "sample_" + type + ".txt";

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

		// export DYLD_LIBRARY_PATH="/opt/intel/mkl/lib/":"/opt/intel/compilers_and_libraries_2016.1.111/mac/compiler/lib":$DYLD_LIBRARY_PATH
		writer.write("export DYLD_LIBRARY_PATH=\"/opt/intel/mkl/lib/\":\"/opt/intel/compilers_and_libraries_2016.1.111/mac/compiler/lib\":$DYLD_LIBRARY_PATH\n");
		writer.write("java -jar -Xmx60g " + this.correlationScript + " filename=" + this.tempName + " writefile=" + covMatFN + " correlation=" + this.correlation + "\n");
		writer.write("cd " + this.writeFolder + "\n");
		// writer.write("pcav0 evd" + covMatFN +"\n");
		writer.write("/Volumes/Promise_RAID/juha/PCA++/pca evd " + covMatFN + "\n");
		writer.write("mv eigenvectors.txt " + this.writeFolder + "GENE.eigenvectors.txt" + "\n");
		writer.write("/Volumes/Promise_RAID/juha/PCA++/pca pc-scores " + type + " " + this.tempName + " GENE.eigenvectors.txt" + "\n");
		writer.write("mv eigenvalues.txt " + this.writeFolder + "GENE.eigenvalues.txt" + "\n");
		writer.write("gzip -f " + covMatFN + "\n");

		// get the number of PCs with Cronbach > 0.7
		String geneEigenVectorFN = "GENE.eigenvectors0.7.txt";
		writer.write("n=0\n" + "while read line; do \n" + "if [[ $line == \"Cronbach\"* ]]; then continue; fi\n" + "compare=$(echo $line'<'0.7 | bc)\n" + "if [[ compare -eq 1 ]]; then break; fi\n" + "((n=$n+1))\n" + "done < cronbach.txt\n" + "echo $n\n" + "cat < GENE.eigenvectors.txt | cut -f1-$n > " + geneEigenVectorFN + "\n");
		writer.write("gzip -f GENE.eigenvectors.txt\n");
		writer.write("gzip -f " + this.tempName + "\n");

		// writer.write("mkdir prediction\n");
		// writer.write(this.javaExe+" "+ this.RNAseqJar + " " +
		// this.getWritePath(getXML("xmlfile")) +"\n");
		// //Create a database from the HPO terms
		//
		// DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd");
		// Date date = new Date();
		// String aucfile =
		// "AUCsAndPValues-"+dateFormat.format(date)+"HPO-unscaled.txt";
		// String hpounscaledfile =
		// "TermGeneZScores-"+dateFormat.format(date)+"HPO-unscaled.txt";
		// String wdir = this.writeFolder+"prediction/";
		// writer.write("node "+ this.populateGenesetDBjs + " " + this.geneDB +
		// " " + this.hpoDB + " " + this.termtype + " " + wdir +
		// this.termDesCoutFile + " " + wdir +"Predictions/"+ hpounscaledfile +
		// " "+ wdir +"Predictions/"+ aucfile + " " + this.writeFolder +
		// "prediction/" + this.geneTermOutFile +"\n");
		// //node
		// /Volumes/Promise_RAID/GeneNetwork/Sipko/05-2016/Tessa/populateGenesetDBTXT.js
		// genedb hpodb HPO
		// /Volumes/Promise_RAID/GeneNetwork/Sipko/06-2016/Juha_QuantNorm_Covariance/Prediction/HPO_terms_22214Samples_Rlog_correlation_+1log.txt
		// /Volumes/Promise_RAID/GeneNetwork/Sipko/06-2016/Juha_QuantNorm_Covariance/Prediction/Predictions/TermGeneZScores-2016-06-10HPO-unscaled.txt
		// /Volumes/Promise_RAID/GeneNetwork/Sipko/06-2016/Juha_QuantNorm_Covariance/Prediction/Predictions/AUCsAndPValues-2016-06-10HPO-unscaled.txt
		// /Volumes/Promise_RAID/GeneNetwork/Sipko/06-2016/Juha_QuantNorm_Covariance/Prediction/HPO_geneterm_22214Samples_Rlog_correlation_+1log.txt
		//
		// //run predictions on 2200 samples
		//
		// writer.write("cd "+ this.tessaPredictionFolder +"\n");
		// writer.write("node " + this.tessaPredictionScript + " " +
		// this.tessaMaxGenesPerTerm + " " + this.getWritePath(this.hpoDB) + " "
		// + this.tessaType + " "+ this.writeFolder+ "\n");
		//
		//
		// //do the Reactiome based prediction
		// //................................
		//
		// //writer.write("java -jar -Xmx60g
		// /Volumes/Promise_RAID/GeneNetwork/Sipko/CorrelationLargeGenes.jar
		// "\n);
		// writeExtraVarsXML(geneEigenVectorFN, date, dateFormat);

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

	private void adjustForPCs(	MatrixStruct inputMatrix,
								int[] PCAadjustments,
								MatrixStruct eigenVectors,
								MatrixStruct scores,
								String writeFolder,
								String vectorFolder,
								MatrixStruct zScores,
								double zScoresCutoff,
								boolean log2,
								boolean correctResultsForSTdevs,
								MatrixStruct chr21,
								int optimalPCremoval,
								String PCs,
								String avgStdevFolder) throws IOException
	{
		// adjustOnZscores();
		MatrixStruct tTestResults = null;
		MatrixStruct difference = null;

		// get PCs to correct for from the user input if any
		ArrayList<Integer> userList = parsePCs(PCs);
		System.out.println("userList=" + userList);
		// JuhaPCA.PCA.log("Calculating variance explained by each PC"); first argument determins the number of principal components to
		// calculate the variance explained for. In this case it is the maximum number of principal components that is corrected for.
		// varianceExplained(PCAadjustments[PCAadjustments.length-1], writeFolder, scores, inputMatrix);//this function only works if the
		// gene standard deviation is set to 1 (and the average of course) JuhaPCA.PCA.log("Calculating variance explained done");

		for (int pcs = -1; pcs < PCAadjustments.length; pcs++)// for all the different numbers of PCs to correct for (e.g. int[] PCAadjustments = new int[]{0,100,2,25,300,500,1000};};)
		{
			MatrixStruct sampleStruct = inputMatrix.copy();// make a copy so you can use the original for when you want to correct for a different number of PCs
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
				if (adjustPCs > eigenVectors.rows())// in case attempting to correct for more PCs than in the eigenvector file, just correct for allt he eigenvectors in the eigenvector files.
					adjustPCs = eigenVectors.rows() + 1;
				writePCName = Integer.toString(adjustPCs);
				for (int p = 1; p < adjustPCs; p++)
				{
					PCsToAdjust.add(p);
				}
			}

			String writeFileName = writeFolder + "PC_1-" + writePCName + "_.txt";
			if (pcs != PCAadjustments.length - 1 || chr21 == null)// if it is not the last batch of PC corrections
			{
				correctPCs(	sampleStruct,
							scores,
							eigenVectors,
							PCsToAdjust,
							null,
							null,
							null,
							optimalPCremoval);
			}
			else // if it is correcting for the full batch of PCs, also calculate how significant the expression changes are for genes on chr 21 for each PC that is subtracted
			{
				correctPCsTrackChr21(	sampleStruct,
										scores,
										eigenVectors,
										PCsToAdjust,
										chr21,
										tTestResults,
										difference,
										optimalPCremoval,
										writeFolder);

				// zScores(writeFileName.replace(".txt", "_Zscores.txt"), sampleStruct);
			}

			sampleStruct = createExtraFiles(sampleStruct,
											correctResultsForSTdevs,
											vectorFolder,
											writeFileName);
			System.out.println("writeFileName =" + writeFileName);
			sampleStruct.write(writeFileName);

		}
	}

	private void correctPCsTrackChr21(	MatrixStruct sampleStruct,
										MatrixStruct scores,
										MatrixStruct eigenVectors,
										ArrayList<Integer> PCsToAdjust,
										MatrixStruct chr21,
										MatrixStruct tTestResults,
										MatrixStruct difference,
										int optimalPCremoval,
										String writeFolder) throws IOException
	{
		if (chr21 == null)
			return;
		chr21.keepRows1Matrix(sampleStruct);

		tTestResults = new MatrixStruct(PCsToAdjust.size(),
										sampleStruct.cols());
		difference = new MatrixStruct(	PCsToAdjust.size(),
										sampleStruct.cols());
		String[] rowHeaders = new String[PCsToAdjust.size()];
		for (int p = 0; p < PCsToAdjust.size(); p++)
		{
			rowHeaders[p] = "PC" + Integer.toString(p + 1);
		}
		tTestResults.setRowHeaders(rowHeaders);
		tTestResults.setColHeaders(sampleStruct.getColHeaders());
		difference.setRowHeaders(rowHeaders);
		difference.setColHeaders(sampleStruct.getColHeaders());

		if (chr21.rows() > 1)
			correctPCs(	sampleStruct,
						scores,
						eigenVectors,
						PCsToAdjust,
						chr21,
						tTestResults,
						difference,
						optimalPCremoval);
		else
		{
			correctPCs(	sampleStruct,
						scores,
						eigenVectors,
						PCsToAdjust,
						null,
						null,
						null,
						optimalPCremoval);
			System.out.println("WARNING!: less then 2 genes of chr21 are present and thus no statistical test can be conducted on the difference after the correction of each PC");
		}

		tTestResults.write(writeFolder + "tTestPerPCchr21.txt");
		difference.write(writeFolder + "differencePerPCchr21.txt");
	}

	private void adjustOnZscores()
	{
		// pca.PCA.log("Copying matrix");
		// MatrixStruct sampleStruct = inputMatrix.copy();
		// pca.PCA.log("Matrix copy done");
		// if(zScoresCutoff >= 0)//if adjusting based on z-scores get the PCs
		// that should be removed
		// {
		// for(int s = 0; s < sampleStruct.cols();s++)
		// {
		// ArrayList<Integer> PCsToAdjust = new ArrayList<Integer>();
		// String outputString = "";
		//
		//
		// correctPCs(sampleStruct, scores, eigenVectors, PCsToAdjust, s);
		//
		//// if(zScores != null)
		//// System.out.println("Sample: "+ s + "/" + zScores.cols()+ "-->" +
		// sampleStruct.getColHeaders()[s] + " Removing " + PCsToAdjust.size() +
		// "/ "+ zScores.rows() +" PCs; Remaining PCs:" + outputString);
		//
		// }
		// String writeFileName = writeFolder+"ZscoreCutoff_" + zScoresCutoff +
		// ".txt";
		// sampleStruct.write(writeFileName);
		// createSmoothedFiles(sampleStruct,correctResultsForSTdevs,vectorFolder,
		// writeFileName);
		// }
	}

	private void varianceExplained(	int i,
									String writeFolder,
									MatrixStruct eigenValues,
									MatrixStruct sampleStruct) throws IOException
	{
		// this only works if the PCscores are calculated based on the matrix where the gene's standard deviation is set to 1

		MatrixStruct explained = new MatrixStruct(	eigenValues.rows(),
													1);
		explained.setRowHeaders(eigenValues.getRowHeaders());
		explained.setColHeaders(sampleStruct.getColHeaders());
		for (int pc = 0; pc < eigenValues.rows(); pc++)// PCs are on the rows
		{
			double PCvar = 0;
			for (int r = 0; r < sampleStruct.rows(); r++)// samples genes are on the rows
			{
				double[] PCscores = eigenValues.getRowValues(pc);
				// System.out.println(PCscores);
				double[] sampleValues = sampleStruct.getRowValues(r);
				// System.out.println(evValues.length + "\t" + sampleValues.length);
				double correlation = Correlation.correlate(	PCscores,
															sampleValues);
				if (pc == 0)
					System.out.println(sampleStruct.getRowHeaders()[r] + "\t" + correlation + "\t" + PCscores[0]);
				double varianceExplained = Math.pow(correlation,
													2);
				PCvar += varianceExplained;
			}
			explained.matrix.set(	pc,
									0,
									PCvar);
		}
		explained.write(writeFolder + "VarianceExplained.txt");
	}

	private void correctPCs(MatrixStruct sampleStruct,
							MatrixStruct scores,
							MatrixStruct eigenVectors,
							ArrayList<Integer> PCsToAdjust,
							MatrixStruct chr21,
							MatrixStruct tTestResults,
							MatrixStruct difference,
							int optimalPCremoval)
	{
		double prevPvalue = 1;
		int outcol = 0;
		for (int pc : PCsToAdjust)// correct this sample for the selected PCs
		{
			if (pc > eigenVectors.rows())
				break;
			for (int s = 0; s < sampleStruct.cols(); s++)
			{
				// remove the signal of this single PC from all the genes
				double[][] out = removeSignalAllgenes(	sampleStruct,
														chr21,
														scores,
														pc,
														s,
														eigenVectors,
														false);

				// check if there is a significant difference between chr21 and the rest
				trackChr21(	out,
							chr21,
							sampleStruct,
							scores,
							pc,
							s,
							eigenVectors,
							optimalPCremoval,
							prevPvalue,
							difference,
							difference,
							tTestResults,
							outcol);
			}
			outcol++;

			if (optimalPCremoval > 0 && pc > 1)// add signal back on that decreases the signal difference between chr21 and the other genes in at least half the down samples
				restoreGoodSignal(	sampleStruct,
									chr21,
									scores,
									pc,
									eigenVectors,
									true,
									optimalPCremoval,
									tTestResults);
		}
	}

	private void restoreGoodSignal(	MatrixStruct sampleStruct,
									MatrixStruct chr21,
									MatrixStruct scores,
									int pc,
									MatrixStruct eigenVectors,
									boolean b,
									int optimalPCremoval,
									MatrixStruct tTestResults)
	{
		// Find out if this last PC should have been kept
		int n = 0;
		for (int c = 0; c < optimalPCremoval; c++)
		{
			if (tTestResults.matrix.get(pc - 1 - 1,
										c) < tTestResults.matrix.get(	pc - 1,
																		c))
				n++;
		}
		if (n >= optimalPCremoval / 2)// the number of samples in which the PC has to be "bad"(decreasing the difference between genes on chr21 and the other genes)
		{
			for (int s = 0; s < sampleStruct.cols(); s++)
			{
				removeSignalAllgenes(	sampleStruct,
										chr21,
										scores,
										pc,
										s,
										eigenVectors,
										true);// add signal back on
			}
		}
	}

	private void trackChr21(double[][] out,
							MatrixStruct chr21,
							MatrixStruct sampleStruct,
							MatrixStruct scores,
							int pc,
							int s,
							MatrixStruct eigenVectors,
							int optimalPCremoval,
							double prevPvalue,
							MatrixStruct difference,
							MatrixStruct difference2,
							MatrixStruct tTestResults,
							int outcol)
	{
		double[] onChr21 = out[0];
		double[] others = out[1];
		if (chr21 != null)
		{
			TTest tTest = new TTest();
			double pValue = tTest.tTest(onChr21,
										others)
					/ 2;
			double avgChr21 = org.apache.commons.math3.stat.StatUtils.mean(onChr21);
			double avgOthers = org.apache.commons.math3.stat.StatUtils.mean(others);
			double diff = avgOthers - avgChr21;

			if (optimalPCremoval == 1)
				if (prevPvalue < pValue)// if the previous p-value is smaller, just add the signal back on
				{
					removeSignalAllgenes(	sampleStruct,
											chr21,
											scores,
											pc,
											s,
											eigenVectors,
											true);// add signal back on
					pValue = prevPvalue;
				}
			prevPvalue = pValue;

			tTestResults.matrix.set(outcol,
									s,
									pValue);
			difference.matrix.set(	outcol,
									s,
									diff);
		}

	}

	private double[][] removeSignalAllgenes(MatrixStruct sampleStruct,
											MatrixStruct chr21,
											MatrixStruct scores,
											int pc,
											int s,
											MatrixStruct eigenVectors,
											boolean add)
	{
		// correct all the genes for this PC
		double[] onChr21 = null, others = null;
		// , beforeOthers = null,Others =null
		if (chr21 != null)
		{
			onChr21 = new double[chr21.rows()];
			others = new double[sampleStruct.rows() - chr21.rows()];
		}
		int chr21Index = 0;
		int othersIndex = 0;
		for (int gene = 0; gene < sampleStruct.rows(); gene++)
		{
			double signal = scores.matrix.get(	pc - 1,
												s)
					* eigenVectors.matrix.get(	pc - 1,
												gene);
			// if(gene == 0 && pc <= 2 && s==0)
			// if(sampleStruct.getRowHeaders()[gene].contains("ENSG00000000003")
			// && s==0)
			// System.out.println("gene = " +sampleStruct.getRowHeaders()[gene]
			// + " geneEigen = " + eigenVectors.getColHeaders()[gene] + " PC ="
			// +pc+" signal = \t" + signal +
			// " score ="+ scores.matrix.get(pc-1, s) + " vectorValue= " +
			// eigenVectors.matrix.get(pc-1, gene)+ "val before = "
			// +sampleStruct.matrix.get(gene, s));
			if (add)// add the signal instead of removing it
				signal *= -1;
			sampleStruct.matrix.add(gene,
									s,
									-signal);
			// if(gene == 0 && pc <= 2 && s==0)
			// if(sampleStruct.getRowHeaders()[gene].contains("ENSG00000000003")
			// && s==0)
			// System.out.println("gene = " +sampleStruct.getRowHeaders()[gene]
			// + " geneEigen = " + eigenVectors.getColHeaders()[gene] + " PC ="
			// +pc+"val after = " +sampleStruct.matrix.get(gene, s));
			if (chr21 != null)
			{
				if (chr21.rowHash.containsKey(sampleStruct.getRowHeaders()[gene]))
				{
					onChr21[chr21Index] = sampleStruct.matrix.get(	gene,
																	s);
					chr21Index++;
				}
				else
				{
					others[othersIndex] = sampleStruct.matrix.get(	gene,
																	s);
					othersIndex++;
				}
			}
		}
		return new double[][] { onChr21, others };
	}

	private MatrixStruct createExtraFiles(	MatrixStruct sampleStruct,
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
		String zscoreWriteFN = new File(writeFileName).getName().replace(	".txt",
																			"_zScores.txt");
		return Zscore.zScores(	this.writeFolderCorrected,
								zscoreWriteFN,
								sampleStruct,
								this.avgStdevFolder,
								zscoreWriteFN.replace(	"_zScores.txt",
														"stats.txt"),
								true);
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

	void keepTopPercentage(	MatrixStruct sampleStruct,
							String averagesFN,
							double topPercent,
							String writeFileName,
							boolean lowest,
							boolean writeAll) throws IOException
	{
		JuhaPCA.PCA.log("20. Highest " + topPercent * 100 + "% only");
		MatrixStruct averages = new MatrixStruct(averagesFN);
		averages.sortCol(0);
		averages.write(averagesFN.replace(	".txt",
											"SAMPLE_Norm_columnAverages_SORTED_ALL.txt"));
		MatrixStruct part = keepPart(	averages,
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

	private void devideBySTdev(	MatrixStruct sampleStruct,
								boolean correctResultsForSTdevs,
								String vectorFolder,
								String writeFileName) throws IOException
	{

		MatrixStruct copy = sampleStruct.copy();
		if (correctResultsForSTdevs)
		{
			JuhaPCA.PCA.log("16. Divide by standard deviation");
			MatrixStruct STdevs = new MatrixStruct(vectorFolder + "gene_STDevs.txt");
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

	private MatrixStruct keepPart(	MatrixStruct averages,
									double topPercent,
									boolean lowest) // returns the last part of the matrix
	{
		int topX = (int) (topPercent * averages.rows());
		MatrixStruct part = new MatrixStruct(	topX,
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
		part.rowHash = MatrixStruct.makeHash(part.getRowHeaders());
		return part;
	}

	private void smoothSignal(	MatrixStruct sampleStruct,
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
		if (PCsToAdjust.contains("null"))
			return null;
		ArrayList<Integer> PCs = new ArrayList<Integer>();
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
		File[] files = new File[7];

		files[0] = new File(this.expFile);
		files[1] = new File(this.geneNameFile);
		files[2] = new File(this.itemSetFile);
		files[3] = new File(this.RNAseqJar);
		files[4] = new File(this.populateGenesetDBjs);
		files[5] = new File(this.tessaPredictionFolder);
		files[6] = new File(this.correlationScript);

		for (File file : files)
		{
			if (file.exists())
				continue;
			System.out.println("THIS FILE/FOLDER DOES NOT EXIST!\n" + file.getAbsolutePath());
			return false;
		}
		return true;
	}

}

// private static MatrixStruct inDirectPCA(MatrixStruct expressionStruct,
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
// MatrixStruct covMat = new MatrixStruct(expressionStruct.getRowHeaders(),
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
// MatrixStruct[] evds = pca.PCA.evd(covMat, Paths.get(writeFolder+"SAMPLE"));
// MatrixStruct eigenVectors = evds[0];
// MatrixStruct PCeigenvalues = evds[1];
//
// pca.PCA.log("23. calculating PCscores over the genes");
// String saveNamePCscoresGene = writeFolder + "GENE_PC.txt";
// MatrixStruct[] PCscoresGenesAndAverages =
// PCA.scores(eigenVectors,expressionStruct, saveNamePCscoresGene);
// MatrixStruct PCscoresGenes = PCscoresGenesAndAverages[0];
// PCscoresGenesAndAverages=null;
// System.gc();
//
// pca.PCA.log("24. Transform PCscores to eigenvectors of Genes");
// String saveNameEigenVectorsOverGenes = writeFolder +
// "GENE.eigenvectors.txt.gz";
// MatrixStruct geneEigenVectors = PCA.transform(PCscoresGenes, PCeigenvalues,
// saveNameEigenVectorsOverGenes);
// return geneEigenVectors;
// }
