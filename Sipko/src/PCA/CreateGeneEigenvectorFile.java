package PCA;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
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
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Text;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import eqtlmappingpipeline.normalization.Normalizer;
import no.uib.cipr.matrix.NotConvergedException;
import pca.MatrixStruct;
import pca.PCA;

public class CreateGeneEigenvectorFile 
{
	/*
	1. load expression matrix
	//2. removeDuplicates (any samples with correlation above <removeDuplicates> threshold)
	3. put genes on rows if they are not already
	4. Sort genes by chromosome and remove any that are not on in the <chromLocationsFile>
	//5. add random values (add a random value of <randomAddValue> to any gene with an expression lower than <randomAddValue>)
	//6. remove low expression (remove genes that do not have an expression of at least <minExpression> in at least <minSamplesExpressed> samples.
	//7. remove a single gene with name <gene> (if you want to)
	//8. quantile normalize 
	//8. Voom normalize
	8. Rlog normalization (without logging data). This is referring to part of the normalization DEseq uses
	9. 2log transform (adds <value> before doing transformation)
	//10. Take the <highestExpressed> % genes
	11. Transpose
	//12. keep genes with highest <TopVariance> % variance
	//13. correct Input for standard deviation if <correctInputForSTdevs> == true
	//14. calculate spearman correlation if <spearman> >=0. Uses <spearman> as minimum expression (if lower, expression is set to this number)
	//15. if <setLowestToAverage> == True, all the lowest values to the average (if allowing log(0), all 0 will become this average))
	16. Calculate column average and store them
	17. Correct for column averages
	18. Calculate rowaverages and store them
	//19. Correct for row averages
	//20. correctInputForSTdevsAfterCenter

	21. Put genes on rows if they are not already
	22. Remove genes that have no variance
	23. calculate correlation/covariance matrix
	24. PCA over matrix
	25. Calculate cronbach alphas (and therefore PCscores over samples)
	Different script(s):
	26. Take all PCs>0.7 and recalculate gene eigenvectors
	27. Put in database 
	*/
	//static CreateGeneEigenvectorFileVariables Var = new CreateGeneEigenvectorFileVariables();
	static Var var = null;
	static Document xml= null;
	
	public static void main(String[] args) throws IOException, NotConvergedException, InterruptedException, ParserConfigurationException, TransformerException
	{
		var = new Var();
		if(args.length == 0)
			checkArgs(args);
		
		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase()){
				case "json":
					var = var.readVars(value);
				break;
				case "filename":
					var.expFile = value;
					break;
				case "chrom":
					var.chromLocationsFile = value;
					break;
				case "writeall":
					var.writeAll = Boolean.parseBoolean(value);
					break;
				case "correctstdevs":
					var.correctInputForSTdevs = Boolean.parseBoolean(value);
					break;
				case "correctstdevsaftercenter":
					var.correctInputForSTdevsAfterCenter = Boolean.parseBoolean(value);
					break;
				case "correcttotalreadcount":
					var.correctTotalReadCount = Double.parseDouble(value);
					break;
				case "rlog":
					var.rLog = Double.parseDouble(value);
					break;
				case "log2":
					var.log2 = Boolean.parseBoolean(value);
					break;
				case "randomvalue":
					var.randomValue = Double.parseDouble(value);
					break;
				case "writefolder":
					var.writeFolder = value;
					break;
				case "duplicate":
					var.duplicateCutoff = Double.parseDouble(value);
					break;
				case "highestexpressed":
					var.highestExpressed = Double.parseDouble(value);
					break;
				case "noqn":
					var.skipQuantileNorm = Boolean.parseBoolean(value);
					break;
				case "stdevcutoff":
					var.STdevCutoff = Boolean.parseBoolean(value);
					break;
				case "zscores":
					var.zScores = Boolean.parseBoolean(value);
					break;
				case "topvariance":
					var.topVariance = Double.parseDouble(value);
					break;
				case "directpca":
					var.directPCA = Boolean.parseBoolean(value);
					break;
				case "spearman":
					var.spearman = Double.parseDouble(value);
					break;
				case "addlogval":
					var.addLogVal = Double.parseDouble(value);
					break;
				case "correlation":
					var.correlation = Boolean.parseBoolean(value);
					break;
				case "lowesttoaverage":
					var.setLowestToAverage = Boolean.parseBoolean(value);
					break;
				case "adjustsampleaverages":
					var.adjustSampleAverages = Boolean.parseBoolean(value);
					break;	
				case "removegene":
					var.removeGene = value;
					break;
				case "minsamplesexpressed":
					var.minSamplesExpressed = Integer.parseInt(value);
					break;
				case "minexpression":
					var.minExpression = Integer.parseInt(value);
					break;
				case "genenamefile":
					var.geneNameFile = value;
					break;
				case "itemsetfile":
					var.itemSetFile = value;
					break;
				case "correctgc":
					var.correctGCSamples = value;
					break;
				default:
					checkArgs(args);
					System.out.println("Incorrect argument supplied; exiting");
					System.exit(1);
			}
		}
		if(var.writeFolder == null)
			var.writeFolder = var.expFile.replace(".txt", "").replace(".gz", "")+"/";
		
		writeParameters();
		var.filePathsExist();
		
		//normalize and center the data
		normalize();
		//1. calculate the correlation/covariance matrix and 
		//2. do the PCA, 
		//3. calculate cronbach alphas
		pca();
	}
	private static void addToXML(String variableName, String variable) 
	{
		addToXML(variableName, variable, false);
	}
	
	private static void addToXML(String variableName, String variable, boolean write) 
	{
		Element config = xml.getDocumentElement();
		Element xmlDate = xml.createElement(variableName.toLowerCase());
		config.appendChild(xmlDate);
		Text stringDate = xml.createTextNode(variable);
			xmlDate.appendChild(stringDate);
		if(write == true)
			writeXML();
	}
	static void checkArgs(String[] args) 
	{
		System.out.println("This script calculates the eigenvectors over the genes and uses the following input:\n"
				+ "fileName=<fileName> - Expression file (samples on rows, gene names on columns)\n"
				+ "chrom=<fileName chromosome file> - File with genes names in 1st column, chromosome in 2nd, position in 3rd (default=null)\n"
				+ "writeAll=<false/true> - writes all files from intermediary steps (default=true)\n"
				+ "correctSTdevs=<false/true> - Corrects for STdevs prior to calculating covariance matrix (default=false)\n"
				+ "log2=<false/true> - Log2 transformation after quantile normalization (default=true)\n"
				+ "randomValue=<number> - adds a random value*<number> to numbers smaller then <number> (default=1)\n"
				+ "writeFolder=<folderName> - folder to write the results in (default=filename.replace(.txt,/)\n"
				+ "duplicate=<number> - correlation cutoff to consider samples as duplicates (default=0.999)\n"
				+ "skipQN=<false/true> - whether to skip the Quantile normalizatin or not (default=false)\n"
				+ "highestExpressed=<number> - Percentage of genes to keep (genes most highly expressed on average)");
		System.exit(1);
	}
	public static void normalize() throws IOException, NotConvergedException, InterruptedException
	{		
		pca.PCA.log(" 1. Reading expression file");
		MatrixStruct expressionStruct = new MatrixStruct(var.expFile);
		
		if(var.writeAll)
			expressionStruct.write(var.writeFolder+"startMatrix.txt.gz");
		
		RemoveDuplicates.removeDuplicates(expressionStruct, var.duplicateCutoff, var.writeFolder, var.writeAll);

		pca.PCA.log(" 3. Transposing");
		expressionStruct.putGenesOnRows();
		
		//sorting on chromosome locations
		expressionStruct = SortChromosome.sort(expressionStruct, var.chromLocationsFile);
		
		//adding random values
		expressionStruct.addRandomValues(var.randomValue);
		expressionStruct.write(var.writeFolder+"randAddedToBelow_" +var.randomValue + ".txt.gz");
		
		//remove genes without variance
		expressionStruct.removeNoVariance(var.writeFolder+"noVarRemoved.txt.gz");
		
		if(var.minExpression > 0)
			expressionStruct.removeLowExpression((var.writeFolder+"RemovedBelow" +var.minExpression+ "inAtLeast" + var.minSamplesExpressed + "samples.txt.gz"), var.minSamplesExpressed, var.minExpression);
		System.out.println("Rows = " + expressionStruct.rows() + " cols = " + expressionStruct.cols());
		
		if(var.removeGene != null)
			expressionStruct.removeRow(var.removeGene);
		
		if(!var.skipQuantileNorm)//quantile Normalize
			QuantileNormalize.quantileNormalize(expressionStruct, var.writeFolder, var.writeAll);
		
		if(var.correctTotalReadCount > 0)
			CorrectReadcounts.correct(var.writeFolder, var.correctTotalReadCount, expressionStruct, var.writeAll, 0.5);
		
		if(var.rLog > 0)//Does not log the values, just does the correction
		{
			String writeGeoFN = var.writeFolder+ "geoMean.txt";
			addToXML("geofile",writeGeoFN);
			writeXML();
			RLog.rLog(var.writeFolder, expressionStruct, var.rLog, var.writeAll, writeGeoFN);
		}
		
		if(var.log2)
			LogTransform.log2(var.writeFolder, expressionStruct, var.writeAll, var.addLogVal);
		
		if(var.highestExpressed >0 && var.highestExpressed < 1)
			HighestExpressed.highestExpressed(expressionStruct, var.skipQuantileNorm, var.correctTotalReadCount, var.writeFolder, var.highestExpressed, var.STdevCutoff, var.writeAll);
		
		pca.PCA.log("11. Transposing");
		expressionStruct.transpose();
		
		if(var.topVariance > 0 && var.topVariance !=1)
			KeepTopVariance.keepTopVariance(expressionStruct, true, var.ignoreLowestValues, var.writeFolder, var.topVariance);
		
		pca.PCA.log("12 Calculating STdevs");
		System.gc();System.gc();
		MatrixStruct stDevs = expressionStruct.stDevCols();
		stDevs.write(var.writeFolder + "gene_STDevs.txt");
		System.out.println("Rows = " + expressionStruct.rows()+ " cols = " + expressionStruct.cols());
		
		if(var.correctInputForSTdevs)
		{
			pca.PCA.log("13 Divide all gene values by STdev for each gene ");	
			expressionStruct.divideBy(stDevs,false);//false corrects columns, true corrects rows
			
			pca.PCA.log("14 Writing matrix divided by gene STdevs");
			expressionStruct.write(var.writeFolder + "_DividedBySGenesSTdev.txt.gz");
		}	
		
		//expressionStruct.isExpressed(null, expressionStruct.cols()*0.05,10);//genes need to have an expression of at least 1 or more in at least 5% of all the samples.
		if(var.spearman >= 0)
		{
			Spearman.ranks(var.writeFolder, expressionStruct, var.spearman);
			var.correlation = true;
		}
		
		if(var.setLowestToAverage)//this will cause the lowest values not to contribute to correlation or covariance
			LowestToAverage.lowestToAverage(expressionStruct);
		
		pca.PCA.log("15. Calculating column averages");
		MatrixStruct colAverages = expressionStruct.getAveragesPerCol();
		String columnavgsfile = var.writeFolder+ "SAMPLE_Norm_columnAverages.txt";
		colAverages.write(columnavgsfile);
		addToXML("columnavgsfile", columnavgsfile);
		writeXML();
		pca.PCA.log("16. Centering: Adjusting for column averages");
		expressionStruct.adjustForAverageAllCols(colAverages);
		expressionStruct.write(var.writeFolder +"SAMPLE_adjustedForGeneAverages.txt.gz");
		pca.PCA.log("17. Calculating row averages");
		MatrixStruct rowAverages = expressionStruct.getAveragesPerRow();
		rowAverages.write(var.writeFolder+ "SAMPLE_Norm_rowAverages.txt");
		if(var.adjustSampleAverages)
		{
			pca.PCA.log("18. Adjusting for row averages");
			expressionStruct.adjustForAverageAllrows(rowAverages);
		}
	
		String expNormLogCentFile = var.writeFolder+"MATRIX_Centered.txt.gz";
		pca.PCA.log("19. Writing centered file in: " + expNormLogCentFile);
		expressionStruct.write(expNormLogCentFile);
		
		double[] cutoffs = KeepTopVariance.getCutoffs(var.ignoreLowestValues, expressionStruct);
		MatrixStruct geneEigenVectors = null;
		
		if(var.correctInputForSTdevsAfterCenter)
		{
			MatrixStruct stDevsRows = expressionStruct.stDevRows();
			stDevsRows.write(var.writeFolder + "_SampleStDevs.txt");
			pca.PCA.log("13 Divide all gene values by STdev for each sample");	
			expressionStruct.divideBy(stDevsRows,true);//false corrects columns, true corrects rows
			
			pca.PCA.log("14 Writing matrix divided by gene STdevs");
			expressionStruct.write(var.writeFolder + "_DividedBySGenesSTdev.txt.gz");
		}
		
		if(var.GCgenes.length()>0)
		{
			MatrixStruct geneGC = new MatrixStruct(var.GCgenes);
			expressionStruct = GCcontent.calculateAndCorrect(expressionStruct, geneGC, var.writeFolder+"gCperSampleWriteFN.txt", var.writeFolder + "GCcorrected.txt.gz");
		}

//		if(var.directPCA)
			directPCA(expressionStruct, var.writeFolder, var.correlation, null);
//		else
//			geneEigenVectors = inDirectPCA(expressionStruct, var.writeFolder, var.correlation, null);
//
//		pca.PCA.log("25. Calculating PCscores for all samples");
//		MatrixStruct[] scoreResults = PCA.scores(geneEigenVectors,expressionStruct, writeFolder+"SAMPLE_PC.txt");
//		MatrixStruct PCsampleScores = scoreResults[0];
//		
//		if(zScores == true)
//		{
//			pca.PCA.log("26. Calculating Z-scores for all PCscores for all samples");
//			MatrixStruct zScoreStats = Zscore.changeToZscores(PCsampleScores);
//			
//			pca.PCA.log("27. Writing Z-scores");
//			PCsampleScores.write(writeFolder+ "pcZscoresSamples.txt");
//			zScoreStats.write(writeFolder+ "pcZscores_Stats.txt");
//		}
//		pca.PCA.log("Files written to: " + writeFolder);
	}
	
//	private static MatrixStruct inDirectPCA(MatrixStruct expressionStruct,
//			String writeFolder, boolean correlation, double[] cutoffs) throws IOException, NotConvergedException {
//		if(expressionStruct.getRowHeaders()[0].contains("ENSG0") && expressionStruct.getRowHeaders()[0].contains("ENST0"))
//		{
//			System.out.println("Genes should be on Cols and are not, transposing (assumes geneIDs start with ENSG0 or ENST0)");
//			expressionStruct.transpose();
//		}
//		String type = "covariance";
//		//expressionStruct.removeNoVariance(writeFolder+"noVarRemoved.txt");
//		
//		if(correlation)
//		{	
//			type = "correlation";
//		}
//		ConcurrentCovariation calculator = new ConcurrentCovariation(20);
//		double[][] inMat = expressionStruct.getMatrix();
//		double[][] covMatrix = calculator.pairwiseCovariation(inMat,false, null, expressionStruct.getRowHeaders(),correlation, cutoffs);
//		
//		MatrixStruct covMat = new MatrixStruct(expressionStruct.getRowHeaders(), expressionStruct.getRowHeaders(), covMatrix);
//		
//		
//		pca.PCA.log("21. Writing covariance matrix over the samples");
//		String covMatFN = writeFolder+"SAMPLE_"+type+".txt.gz";
//		covMat.write(covMatFN);
//		
//		//to test something
//		Matrix covMatOtherFormat = new Matrix(expressionStruct.getRowHeaders(), expressionStruct.getRowHeaders(), covMatrix);
//		covMatOtherFormat.getAverageCols(true)
//		  .write(covMatFN.replace(".txt", "_Absolute_averages.txt"));
//		covMatOtherFormat.write(covMatFN);
//		
//		inMat = null; covMatrix = null; System.gc();System.gc();
//				
//		pca.PCA.log("22. calculating eigenvalues over the samples");
//		MatrixStruct[] evds = pca.PCA.evd(covMat, Paths.get(writeFolder+"SAMPLE"));
//		MatrixStruct eigenVectors = evds[0];
//		MatrixStruct PCeigenvalues = evds[1];
//		
//		pca.PCA.log("23. calculating PCscores over the genes");
//		String saveNamePCscoresGene = writeFolder + "GENE_PC.txt";
//		MatrixStruct[] PCscoresGenesAndAverages = PCA.scores(eigenVectors,expressionStruct, saveNamePCscoresGene);
//		MatrixStruct PCscoresGenes = PCscoresGenesAndAverages[0];
//		PCscoresGenesAndAverages=null;
//		System.gc();
//		
//		pca.PCA.log("24. Transform PCscores to eigenvectors of Genes");
//		String saveNameEigenVectorsOverGenes = writeFolder + "GENE.eigenvectors.txt.gz";
//		MatrixStruct geneEigenVectors = PCA.transform(PCscoresGenes, PCeigenvalues, saveNameEigenVectorsOverGenes);
//		return geneEigenVectors;
//	}
	private static void directPCA(MatrixStruct expressionStruct, String writeFolder, boolean correlation, double[] cutoffs) throws IOException, InterruptedException {
		expressionStruct.putGenesOnRows();
		
		ConcurrentCovariation calculatorGenes = new ConcurrentCovariation(20);
		String type = "covariance";
		//expressionStruct.removeNoVariance(writeFolder+"noVarRemoved.txt");
		
		if(correlation)	
			type = "correlation";
		
		addToXML("type", type);
		String tempName = writeFolder+"pre_Correlation_Or_Covariance.txt";
		expressionStruct.write(tempName);
		addToXML("normalizedExpressionFile", tempName);
		expressionStruct = null; System.gc();System.gc();
		String coExpMatFN = writeFolder+"gene_"+type+".txt";
		addToXML("coexpmatfile", coExpMatFN, true);
	}
	private static void writeParameters() throws IOException, ParserConfigurationException, TransformerException {
		makeFolder(var.writeFolder);
		
		DateFormat dateFormat = new SimpleDateFormat("dd/MM/yyyy HH:mm:ss");
		Date date = new Date();
		//this can be removed since it makes an .xml file now
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(var.writeFolder+"parameters.txt")));
		
		writer.write("Date\t" + dateFormat.format(date) + "\n");
		writer.write("Input file\t" + var.expFile + "\n");
		writer.write("writeFolder\t" + var.writeFolder + "\n");
		writer.write("chromLocationsFile\t" + var.chromLocationsFile + "\n");
		writer.write("writeAll\t" + var.writeAll + "\n");
		writer.write("correctTotalReadCount\t" + var.correctTotalReadCount + "\n");
		writer.write("log2\t" + var.log2 + "\n");
		writer.write("correctInputForSTdevs\t" + var.correctInputForSTdevs + "\n");
		writer.write("randomValue\t" + var.randomValue + "\n");
		writer.write("duplicateCutoff\t" + var.duplicateCutoff + "\n");
		writer.write("TopXhighestExpressed\t" + var.highestExpressed + "\n");
		writer.write("skipQuantileNorm\t" + var.skipQuantileNorm + "\n");
		writer.write("STdevCutoff\t" + var.STdevCutoff + "\n");
		writer.write("zScores\t" + var.zScores + "\n");
		writer.write("directPCA\t" + var.directPCA + "\n");
		writer.write("spearman\t" + var.spearman + "\n");
		writer.write("setLowestToAverage\t" + var.setLowestToAverage + "\n");
		writer.write("adjustSampleAverages\t" + var.adjustSampleAverages + "\n");
		writer.write("removeGene\t" + var.removeGene + "\n");
		writer.write("minExpression\t" + var.minExpression + "\n");
		writer.write("minSamplesExpressed\t" + var.minSamplesExpressed + "\n");
		writer.write("correlation\t" + var.correlation + "\n");
		writer.write("addLogVal\t" + var.addLogVal + "\n");
		
		writer.close();
		//write an XML file:
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		DocumentBuilder db = dbf.newDocumentBuilder();
		xml = db.newDocument();
			
		Element config = xml.createElement("config");
		xml.appendChild(config);
		addToXML("date", dateFormat.format(date));
		addToXML("expressionfile", var.expFile);
		
		addToXML("writeFolder", var.writeFolder);
		addToXML("chromLocationsFile", var.chromLocationsFile);
		addToXML("writeAll", Boolean.toString(var.writeAll));
		addToXML("correctTotalReadCount", Double.toString(var.correctTotalReadCount));
		addToXML("log2", Boolean.toString(var.log2));

		addToXML("correctInputForSTdevs", Boolean.toString(var.correctInputForSTdevs));
		addToXML("randomValue", Double.toString(var.randomValue));
		addToXML("duplicateCutoff", Double.toString(var.duplicateCutoff));
		addToXML("TopXhighestExpressed", Double.toString(var.highestExpressed));
		addToXML("skipQuantileNorm", Boolean.toString(var.skipQuantileNorm));
		addToXML("STdevCutoff", Boolean.toString(var.STdevCutoff));
		addToXML("zScores", Boolean.toString(var.zScores));
		addToXML("directPCA", Boolean.toString(var.directPCA));

		addToXML("spearman", Double.toString(var.spearman));
		addToXML("setLowestToAverage", Boolean.toString(var.setLowestToAverage));
		addToXML("adjustSampleAverages", Boolean.toString(var.adjustSampleAverages));
//		addElement("removeGene", removeGene);
		addToXML("minExpression", Integer.toString(var.minExpression));
		addToXML("minSamplesExpressed", Integer.toString(var.minSamplesExpressed));
		addToXML("correlation", Boolean.toString(var.correlation));
		addToXML("addLogVal", Double.toString(var.addLogVal));
		addToXML("xmlfile", var.xmlFN);
		addToXML("correctGC", var.correctGCSamples);
		writeXML();
		var.jsonFN=var.writeFolder+"config.json";
		if(var.jsonFN == null)
			System.exit(1);
		var.writeVars(var.jsonFN);
	}
	
	private static void writeXML() {
		String writeFN = var.writeFolder+getXML("xmlfile");
		TransformerFactory tff = TransformerFactory.newInstance();
		Transformer tf;
		try {
			tf = tff.newTransformer();
		
			DOMSource source = new DOMSource(xml);
			StreamResult result = new StreamResult(new File(writeFN).getAbsolutePath());
			tf.setOutputProperty(OutputKeys.INDENT, "yes");
			tf.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "2");
			tf.transform(source, result);
		} catch (TransformerConfigurationException e) {	e.printStackTrace();} catch (TransformerException e) {e.printStackTrace();}
	}

	static void makeFolder(String writeFolder) 
	{
		File folder = new File(writeFolder);
		if(!folder.exists())
		{
			folder.mkdir();
		}
		
	}
	static void pca() throws IOException, InterruptedException
	{	
		System.gc();System.gc();System.gc();System.gc();
		String tempName = getXML("normalizedexpressionfile");
		String covMatFN = getXML("coexpmatfile");
		String writeFolder = getXML("writefolder");
		String type = getXML("type");
		
		Runtime runtime = Runtime.getRuntime();
		NumberFormat format = NumberFormat.getInstance();
		long maxMemory = runtime.maxMemory();
		long allocatedMemory = runtime.totalMemory();
		long freeMemory = runtime.freeMemory();

		System.out.println("free memory: " + format.format(freeMemory / 1024/1024/1024) + "<br/>");
		System.out.println("allocated memory: " + format.format(allocatedMemory / 1024/1024/1024) + "<br/>");
		System.out.println("max memory: " + format.format(maxMemory / 1024/1024/1024) + "<br/>");
		System.out.println("total free memory: " + format.format((freeMemory + (maxMemory - allocatedMemory)/1024/1024/1024)));
		
		PCA.log("Matrix decomposition");
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(writeFolder +"export.sh")));
		
		addToXML("correlationScript",var.correlationScript);
		
		//export DYLD_LIBRARY_PATH="/opt/intel/mkl/lib/":"/opt/intel/compilers_and_libraries_2016.1.111/mac/compiler/lib":$DYLD_LIBRARY_PATH
		writer.write("export DYLD_LIBRARY_PATH=\"/opt/intel/mkl/lib/\":\"/opt/intel/compilers_and_libraries_2016.1.111/mac/compiler/lib\":$DYLD_LIBRARY_PATH\n");
		writer.write("java -jar -Xmx60g "+ var.correlationScript +" filename=" + tempName + " writefile=" + covMatFN + " correlation=" + var.correlation +"\n");
		writer.write("cd " + writeFolder +"\n");
		writer.write("/Volumes/Promise_RAID/juha/PCA++/pca evd " + covMatFN +"\n");
		writer.write("mv eigenvectors.txt " + writeFolder + "GENE.eigenvectors.txt"+"\n");
		writer.write("/Volumes/Promise_RAID/juha/PCA++/pca pc-scores " + type + " " + tempName + " GENE.eigenvectors.txt" + "\n");
		writer.write("mv eigenvalues.txt " + writeFolder + "GENE.eigenvalues.txt"+"\n");
		writer.write("gzip " + covMatFN + "\n");
		
		//get the number of PCs with Cronbach > 0.7
		String geneEigenVectorFN = "GENE.eigenvectors0.7.txt";
		writer.write("n=0\n"
				+ "while read line; do \n"
				+ "if [[ $line == \"Cronbach\"* ]]; then continue; fi\n"
				+ "compare=$(echo $line'<'0.7 | bc)\n"
				+ "if [[ compare -eq 1 ]]; then break; fi\n"
				+ "((n=$n+1))\n"
				+ "done < cronbach.txt\n"
				+ "echo $n\n"
				+ "cat < GENE.eigenvectors.txt | cut -f1-$n > "+geneEigenVectorFN+"\n");
		writer.write("gzip GENE.eigenvectors.txt\n");
		writer.write("gzip " + tempName + "\n");
		
		//do the HPO term based prediction
		addToXML("wdir",var.writeFolder+"prediction/");
		addToXML("datafile",var.getWritePath(geneEigenVectorFN));
		addToXML("genenamefile",var.geneNameFile);
		addToXML("itemsetfile",var.itemSetFile);
		addToXML("genetermoutfile", var.geneTermOutFile);
		addToXML("termdescoutfile", var.termDesCoutFile);
		addToXML("itemtype",var.itemType);
		addToXML("limittoitemsingenesetfile",var.limitToItemsInGenesetFile);
		addToXML("maxitems",var.maxItems);
		addToXML("minitems",var.minItems);
		addToXML("absminitems",var.absMinItems);
		addToXML("usettest",var.useTtest);
		addToXML("runrealanalysis",var.runRealAnalysis);		
		addToXML("nrpermutations",var.nrPermutations);
		addToXML("writebinary",var.writeBinary);
		addToXML("logtofile",var.logToFile);
		addToXML("label",var.label);
		addToXML("javaexe",var.javaExe);
		addToXML("rnaseqjar",var.RNAseqJar);
		
		
		writer.write("mkdir prediction\n");
		writer.write(var.javaExe+" "+ var.RNAseqJar + " " + var.getWritePath(getXML("xmlfile")) +"\n");	
		//Create a database from the HPO terms
		
		addToXML("populategenesetdbjs",var.populateGenesetDBjs);
		addToXML("genedb",var.geneDB);
		addToXML("hpodb",var.hpoDB);
		addToXML("termtype",var.termtype);
		DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd");
		Date date = new Date();		
		addToXML("hpounscaledfile", "TermGeneZScores-"+dateFormat.format(date)+"HPO-unscaled.txt");
		addToXML("aucfile", "AUCsAndPValues-"+dateFormat.format(date)+"HPO-unscaled.txt");
		
		writer.write("node "+ var.populateGenesetDBjs + " " + getXML("genedb") + " " + getXML("hpodb") + " " + getXML("termtype") + " " + getXML("wdir")+getXML("termdescoutfile")+ " "+ getXML("wdir")+"Predictions/"+getXML("hpounscaledfile") + " "+ getXML("wdir")+"Predictions/"+getXML("aucfile") + " " + getXML("wdir")+getXML("genetermoutfile") +"\n");
		//node /Volumes/Promise_RAID/GeneNetwork/Sipko/05-2016/Tessa/populateGenesetDBTXT.js genedb hpodb HPO /Volumes/Promise_RAID/GeneNetwork/Sipko/06-2016/Juha_QuantNorm_Covariance/Prediction/HPO_terms_22214Samples_Rlog_correlation_+1log.txt /Volumes/Promise_RAID/GeneNetwork/Sipko/06-2016/Juha_QuantNorm_Covariance/Prediction/Predictions/TermGeneZScores-2016-06-10HPO-unscaled.txt /Volumes/Promise_RAID/GeneNetwork/Sipko/06-2016/Juha_QuantNorm_Covariance/Prediction/Predictions/AUCsAndPValues-2016-06-10HPO-unscaled.txt /Volumes/Promise_RAID/GeneNetwork/Sipko/06-2016/Juha_QuantNorm_Covariance/Prediction/HPO_geneterm_22214Samples_Rlog_correlation_+1log.txt
		
		//run predictions on 2200 samples
		addToXML("tessapredictionfolder",var.tessaPredictionFolder);
		addToXML("tessapredictionscript",var.tessaPredictionScript);
		addToXML("tessamaxgenesperterm",var.tessaMaxGenesPerTerm);
		addToXML("tessatype",var.tessaType);
		writer.write("cd "+ getXML("tessapredictionfolder") +"\n");
		writer.write("node " + getXML("tessapredictionscript") + " " + getXML("tessamaxgenesperterm") + " " + var.getWritePath(getXML("hpodb")) + " " + getXML("tessatype") + " "+ var.writeFolder+ "\n");
		writeXML();
		
		//do the Reactiome based prediction
		//................................
		
		//writer.write("java -jar -Xmx60g /Volumes/Promise_RAID/GeneNetwork/Sipko/CorrelationLargeGenes.jar "\n);
		writer.close();		
		String command = "sh "+ writeFolder + "export.sh";
		run(command);
//		Process p;
//		System.out.println("Starting");
//		p = Runtime.getRuntime().exec(command);
//		p.waitFor();
//		System.out.println("Finished");
		
		if(new File(var.sampleFile).exists() && var.PCs.length()>0)
		{
			var.writeVars(var.jsonFN);
			PCcorrection.main(new String[]{"json="+var.jsonFN});
		}
		System.out.println("All done!");
	}
	private static void run(String command) {
		try
        {            
            Runtime rt = Runtime.getRuntime();
            Process proc = rt.exec(command);
            InputStream stderr = proc.getErrorStream();
            InputStreamReader isr = new InputStreamReader(stderr);
            BufferedReader br = new BufferedReader(isr);
            String line = null;
            while ( (line = br.readLine()) != null)
                System.out.println(line);
            int exitVal = proc.waitFor();
            System.out.println("Process exitValue: " + exitVal);
        } catch (Throwable t)
        {
        	t.printStackTrace();
        }
	}
	private static String getXML(String string) {
		return xml.getElementsByTagName(string).item(0).getTextContent();
	}
}
