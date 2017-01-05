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
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Text;

import JuhaPCA.PCA;
import Tools.ExecCommand;
import eqtlmappingpipeline.normalization.Normalizer;
import no.uib.cipr.matrix.NotConvergedException;

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
	public static Var var = null;
	static Document xml= null;
	public static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko") && args.length < 1)
			return;
		if(args.length == 0)
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
					+ "skipQN=<false/true> - whether to skip the Quantile normalizatin or not (default=true)\n"
					+ "DESeqNorm=<false/true> - whether to use DESeq normalization (default=true)\n"
					+ "highestExpressed=<number> - Percentage of genes to keep (genes most highly expressed on average)"
					+ "genestoinclude=<filename> - File with in the first column the genes to keep (assumes column headers)");
			System.exit(1);
		}
		
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
				case "deseqnorm":
					var.rLog = Boolean.parseBoolean(value);
					break;
				case "rlog":
					var.rLog = Boolean.parseBoolean(value);
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
					var.centerSamples = Boolean.parseBoolean(value);
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
					var.GCgenes = value;
					break;
				case "genestoinclude":
					var.genesToInclude = value;
					break;
				default:
					System.out.println("Incorrect argument supplied:\n"+ args[a] +"\nexiting");
					System.exit(1);
			}
		}
	}
	public static void main(String[] args) throws IOException, NotConvergedException, InterruptedException, ParserConfigurationException, TransformerException
	{
		var = new Var();
		checkArgs(args);
		
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
	private static void selectGenes(MatrixStruct expressionStruct, String genesToInclude) throws IOException {
		if(genesToInclude!=null)
			new MatrixStruct(genesToInclude).keepRows(expressionStruct); //keeps them in the original order
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
	
	public static void normalize() throws IOException, NotConvergedException, InterruptedException
	{		
		JuhaPCA.PCA.log(" 1. Reading expression file");
		MatrixStruct expressionStruct = new MatrixStruct(var.expFile);
		//transposes matrix if genes/transcripts are not on rows
		expressionStruct.putGenesOnRows();
				
		//keep only a subset of genes
		selectGenes(expressionStruct, var.genesToInclude);
		
		if(var.writeAll)
			expressionStruct.write(var.writeFolder+"startMatrix.txt.gz");
				
		//removes duplicates if var.duplicateCutoff <1
		RemoveDuplicates.removeDuplicates(expressionStruct, var.duplicateCutoff, var.writeFolder, var.writeAll);
	
		//sorting on chromosome locations
		expressionStruct = SortChromosome.sort(expressionStruct, var.chromLocationsFile);
		
		//adding random values, if var.randomValue =< 0, skips  this
		expressionStruct.addRandomValues(var.randomValue, var);
		
		System.out.println("Rows1 = " + expressionStruct.rows()+ " cols = " + expressionStruct.cols());
		//remove genes without variance
		
		expressionStruct.removeNoVariance(var.writeFolder+"noVarRemoved.txt.gz");
		//remove genes that have an expression below var.minExpression
		if(var.minExpression > 0)//if var.minExpression is zero or smaller, skip this
			expressionStruct.removeLowExpression((var.writeFolder+"RemovedBelow" +var.minExpression+ "inAtLeast" + var.minSamplesExpressed + "samples.txt.gz"), var.minSamplesExpressed, var.minExpression);
		System.out.println("Rows2 = " + expressionStruct.rows()+ " cols = " + expressionStruct.cols());
		if(var.removeGene != null)//remove the gene defined by var.removeGene
			expressionStruct.removeRow(var.removeGene);
		
		if(!var.skipQuantileNorm)//quantile Normalize if false
			QuantileNormalize.quantileNormalize(expressionStruct, var.writeFolder, var.writeAll);
		
		if(var.correctTotalReadCount > 0)//correct for total number of reads
			CorrectReadcounts.correct(var.writeFolder, var.correctTotalReadCount, expressionStruct, var.writeAll, 0.5);
		
		if(var.rLog)//Does not log the values, just does the DEseq based correction
		{
			//String writeGeoFN = var.writeFolder+ "geoMean.txt";
			RLog.rLog(var.writeFolder, expressionStruct, var.writeAll, null);
		}
		
		if(var.log2)//log2 the data and add <var.addLogVal> before
			LogTransform.log2(var.writeFolder, expressionStruct, var.writeAll, var.addLogVal);
		
		if(var.highestExpressed >0 && var.highestExpressed < 1)//keep the highest <var.highestExpressed> genes only
			HighestExpressed.highestExpressed(expressionStruct, var.skipQuantileNorm, var.correctTotalReadCount, var.writeFolder, var.highestExpressed, var.STdevCutoff, var.writeAll);
		
		JuhaPCA.PCA.log("12 Calculating STdevs");
		System.gc();System.gc();
		MatrixStruct stDevs = expressionStruct.stDevRows();
		stDevs.write(var.writeFolder + "gene_STDevs.txt");
		System.out.println("Rows = " + expressionStruct.rows()+ " cols = " + expressionStruct.cols());
		
		if(var.correctInputForSTdevs)
		{
			JuhaPCA.PCA.log("13 Divide all gene values by STdev for each gene ");	
			expressionStruct.divideBy(stDevs,false);//false corrects columns, true corrects rows
			
			JuhaPCA.PCA.log("14 Writing matrix divided by gene STdevs");
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
		if(var.GCgenes != null)
		{
			MatrixStruct geneGC = new MatrixStruct(var.GCgenes);
			expressionStruct = GCcontent.calculateAndCorrect(expressionStruct, geneGC, var.writeFolder+"gCperSampleWriteFN.txt", var.writeFolder + "GCcorrected.txt.gz");
		}
		
		//expressionStruct.removeNoVariance(var.writeFolder+"noVarRemoved2.txt.gz");
		
		JuhaPCA.PCA.log("15. Calculating column averages");
		MatrixStruct colAverages = expressionStruct.getAveragesPerCol();
		String columnavgsfile = var.writeFolder+ "SAMPLE_Norm_sampleAverages.txt";
		colAverages.write(columnavgsfile);
		if(var.centerSamples)
		{
			JuhaPCA.PCA.log("16. Centering samples: Adjusting for sample averages");
			expressionStruct.adjustForAverageAllSamples(colAverages);
			expressionStruct.write(var.writeFolder +"SAMPLE_adjustedForSampleAverages.txt.gz");
		}
		JuhaPCA.PCA.log("17. Calculating row averages");
		MatrixStruct rowAverages = expressionStruct.getAveragesPerRow();
		rowAverages.write(var.writeFolder+ "SAMPLE_Norm_GeneAverages.txt");
		if(var.centerGenes)
		{
			JuhaPCA.PCA.log("18. Centegering genes: Adjusting for gene averages");
			expressionStruct.adjustForAverageAllGenes(rowAverages);
		}
	
		String expNormLogCentFile = var.writeFolder+"MATRIX_Centered.txt.gz";
		JuhaPCA.PCA.log("19. Writing centered file in: " + expNormLogCentFile);
		expressionStruct.write(expNormLogCentFile);
		
		if(var.correctInputForSTdevsAfterCenter)
		{
			MatrixStruct stDevsRows = expressionStruct.stDevRows();
			stDevsRows.write(var.writeFolder + "_SampleStDevs.txt");
			JuhaPCA.PCA.log("13 Divide all gene values by STdev for each sample");	
			expressionStruct.divideBy(stDevsRows,true);//false corrects columns, true corrects rows
			
			JuhaPCA.PCA.log("14 Writing matrix divided by gene STdevs");
			expressionStruct.write(var.writeFolder + "_DividedBySGenesSTdev.txt.gz");
		}
			
		var.tempName = var.writeFolder+"pre_Correlation_Or_Covariance.txt";
		expressionStruct.write(var.tempName);
		expressionStruct = null; System.gc();System.gc();System.gc();
	}
	
	public static void writeParameters() throws IOException, ParserConfigurationException, TransformerException {
		makeFolder(var.writeFolder);
		
		DateFormat dateFormat = new SimpleDateFormat("dd/MM/yyyy HH:mm:ss");
		Date date = new Date();
		//create an XML file: //need to fix this so I can use it for calling Tessa's/Juha's prediction program that calculates AUCs per Reactome/GO term
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		DocumentBuilder db = dbf.newDocumentBuilder();
		xml = db.newDocument();
		Element config = xml.createElement("config");
		xml.appendChild(config);
		addToXML("date", dateFormat.format(date));
		addToXML("xmlfile", var.xmlFN);
		writeXML();
		//create json File
		var.jsonFN=var.writeFolder+"config.json";
//		if(var.jsonFN == null)
//			System.exit(1);
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
		
		String type = "covariance";
		if(var.correlation)	
			type = "correlation";
		String covMatFN = var.writeFolder+"gene_"+type+".txt";
		
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
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(var.writeFolder +"export.sh")));
	
		//export DYLD_LIBRARY_PATH="/opt/intel/mkl/lib/":"/opt/intel/compilers_and_libraries_2016.1.111/mac/compiler/lib":$DYLD_LIBRARY_PATH
		writer.write("export DYLD_LIBRARY_PATH=\"/opt/intel/mkl/lib/\":\"/opt/intel/compilers_and_libraries_2016.1.111/mac/compiler/lib\":$DYLD_LIBRARY_PATH\n");
		writer.write("java -jar -Xmx60g "+ var.correlationScript +" filename=" + var.tempName + " writefile=" + covMatFN + " correlation=" + var.correlation +"\n");
		writer.write("cd " + var.writeFolder +"\n");
		//writer.write("pcav0 evd" + covMatFN +"\n");
		writer.write("/Volumes/Promise_RAID/juha/PCA++/pca evd " + covMatFN +"\n");
		writer.write("mv eigenvectors.txt " + var.writeFolder + "GENE.eigenvectors.txt"+"\n");
		writer.write("/Volumes/Promise_RAID/juha/PCA++/pca pc-scores " + type + " " + var.tempName + " GENE.eigenvectors.txt" + "\n");
		writer.write("mv eigenvalues.txt " + var.writeFolder + "GENE.eigenvalues.txt"+"\n");
		writer.write("gzip -f " + covMatFN + "\n");
		
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
		writer.write("gzip -f GENE.eigenvectors.txt\n");
		writer.write("gzip -f " + var.tempName + "\n");
		
//		writer.write("mkdir prediction\n");
//		writer.write(var.javaExe+" "+ var.RNAseqJar + " " + var.getWritePath(getXML("xmlfile")) +"\n");	
//		//Create a database from the HPO terms
//		
//		DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd");
//		Date date = new Date();	
//		String aucfile = "AUCsAndPValues-"+dateFormat.format(date)+"HPO-unscaled.txt";
//		String hpounscaledfile = "TermGeneZScores-"+dateFormat.format(date)+"HPO-unscaled.txt";
//		String wdir = var.writeFolder+"prediction/";
//		writer.write("node "+ var.populateGenesetDBjs + " " + var.geneDB + " " + var.hpoDB + " " + var.termtype + " " + wdir + var.termDesCoutFile + " " + wdir +"Predictions/"+ hpounscaledfile + " "+ wdir +"Predictions/"+ aucfile + " " + var.writeFolder + "prediction/" + var.geneTermOutFile +"\n");
//		//node /Volumes/Promise_RAID/GeneNetwork/Sipko/05-2016/Tessa/populateGenesetDBTXT.js genedb hpodb HPO /Volumes/Promise_RAID/GeneNetwork/Sipko/06-2016/Juha_QuantNorm_Covariance/Prediction/HPO_terms_22214Samples_Rlog_correlation_+1log.txt /Volumes/Promise_RAID/GeneNetwork/Sipko/06-2016/Juha_QuantNorm_Covariance/Prediction/Predictions/TermGeneZScores-2016-06-10HPO-unscaled.txt /Volumes/Promise_RAID/GeneNetwork/Sipko/06-2016/Juha_QuantNorm_Covariance/Prediction/Predictions/AUCsAndPValues-2016-06-10HPO-unscaled.txt /Volumes/Promise_RAID/GeneNetwork/Sipko/06-2016/Juha_QuantNorm_Covariance/Prediction/HPO_geneterm_22214Samples_Rlog_correlation_+1log.txt
//		
//		//run predictions on 2200 samples
//		
//		writer.write("cd "+ var.tessaPredictionFolder +"\n");
//		writer.write("node " + var.tessaPredictionScript + " " + var.tessaMaxGenesPerTerm + " " + var.getWritePath(var.hpoDB) + " " + var.tessaType + " "+ var.writeFolder+ "\n");
//		
//		
//		//do the Reactiome based prediction
//		//................................
//		
//		//writer.write("java -jar -Xmx60g /Volumes/Promise_RAID/GeneNetwork/Sipko/CorrelationLargeGenes.jar "\n);
//		writeExtraVarsXML(geneEigenVectorFN, date, dateFormat);
		
		writer.close();	
		String command = "sh "+ var.writeFolder + "export.sh";
		run(command);
		
		if(new File(var.sampleFile).exists() && var.PCs.length()>0)
		{
			var.writeVars(var.jsonFN);
			PCcorrection.main(new String[]{"json="+var.jsonFN});
		}

		System.out.println("All done!");
	}
	private static void writeExtraVarsXML(String geneEigenVectorFN, Date date, DateFormat dateFormat) {
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
		addToXML("adjustSampleAverages", Boolean.toString(var.centerSamples));
//		addElement("removeGene", removeGene);
		addToXML("minExpression", Integer.toString(var.minExpression));
		addToXML("minSamplesExpressed", Integer.toString(var.minSamplesExpressed));
		addToXML("correlation", Boolean.toString(var.correlation));
		addToXML("addLogVal", Double.toString(var.addLogVal));
		
		addToXML("correlationScript",var.correlationScript);
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
		addToXML("populategenesetdbjs",var.populateGenesetDBjs);
		addToXML("genedb",var.geneDB);
		addToXML("hpodb",var.hpoDB);
		addToXML("termtype",var.termtype);
			
		addToXML("hpounscaledfile", "TermGeneZScores-"+dateFormat.format(date)+"HPO-unscaled.txt");
		addToXML("aucfile", "AUCsAndPValues-"+dateFormat.format(date)+"HPO-unscaled.txt");
		addToXML("tessapredictionfolder",var.tessaPredictionFolder);
		addToXML("tessapredictionscript",var.tessaPredictionScript);
		addToXML("tessamaxgenesperterm",var.tessaMaxGenesPerTerm);
		addToXML("tessatype",var.tessaType);
		writeXML();
	}
	private static void run(String command) {
		System.out.println("Shellcommand = " + command);
		ExecCommand exec = new ExecCommand(command);
		System.out.println("execute output: \n"+exec.getOutput());
		System.out.println("execute error: \n"+exec.getError());
	}
	private static String getXML(String string) {
		return xml.getElementsByTagName(string).item(0).getTextContent();
	}
}







//private static MatrixStruct inDirectPCA(MatrixStruct expressionStruct,
//String writeFolder, boolean correlation, double[] cutoffs) throws IOException, NotConvergedException {
//if(expressionStruct.getRowHeaders()[0].contains("ENSG0") && expressionStruct.getRowHeaders()[0].contains("ENST0"))
//{
//System.out.println("Genes should be on Cols and are not, transposing (assumes geneIDs start with ENSG0 or ENST0)");
//expressionStruct.transpose();
//}
//String type = "covariance";
////expressionStruct.removeNoVariance(writeFolder+"noVarRemoved.txt");
//
//if(correlation)
//{	
//type = "correlation";
//}
//ConcurrentCovariation calculator = new ConcurrentCovariation(20);
//double[][] inMat = expressionStruct.getMatrix();
//double[][] covMatrix = calculator.pairwiseCovariation(inMat,false, null, expressionStruct.getRowHeaders(),correlation, cutoffs);
//
//MatrixStruct covMat = new MatrixStruct(expressionStruct.getRowHeaders(), expressionStruct.getRowHeaders(), covMatrix);
//
//
//pca.PCA.log("21. Writing covariance matrix over the samples");
//String covMatFN = writeFolder+"SAMPLE_"+type+".txt.gz";
//covMat.write(covMatFN);
//
////to test something
//Matrix covMatOtherFormat = new Matrix(expressionStruct.getRowHeaders(), expressionStruct.getRowHeaders(), covMatrix);
//covMatOtherFormat.getAverageCols(true)
//.write(covMatFN.replace(".txt", "_Absolute_averages.txt"));
//covMatOtherFormat.write(covMatFN);
//
//inMat = null; covMatrix = null; System.gc();System.gc();
//	
//pca.PCA.log("22. calculating eigenvalues over the samples");
//MatrixStruct[] evds = pca.PCA.evd(covMat, Paths.get(writeFolder+"SAMPLE"));
//MatrixStruct eigenVectors = evds[0];
//MatrixStruct PCeigenvalues = evds[1];
//
//pca.PCA.log("23. calculating PCscores over the genes");
//String saveNamePCscoresGene = writeFolder + "GENE_PC.txt";
//MatrixStruct[] PCscoresGenesAndAverages = PCA.scores(eigenVectors,expressionStruct, saveNamePCscoresGene);
//MatrixStruct PCscoresGenes = PCscoresGenesAndAverages[0];
//PCscoresGenesAndAverages=null;
//System.gc();
//
//pca.PCA.log("24. Transform PCscores to eigenvectors of Genes");
//String saveNameEigenVectorsOverGenes = writeFolder + "GENE.eigenvectors.txt.gz";
//MatrixStruct geneEigenVectors = PCA.transform(PCscoresGenes, PCeigenvalues, saveNameEigenVectorsOverGenes);
//return geneEigenVectors;
//}
