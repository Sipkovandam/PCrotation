package PCA;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.concurrent.ConcurrentCorrelation;
import umcg.genetica.math.stats.concurrent.ConcurrentCovariation;
import java.io.IOException;
import java.nio.file.Paths;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.Hashtable;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

import no.uib.cipr.matrix.NotConvergedException;
import pca.MatrixStruct;
import pca.PCA;

public class CreateGeneEigenvectorFile 
{
	static boolean setLowestToAverage = false;
	static boolean adjustSampleAverages = true;
	static String removeGene = null;
	static int minSamplesExpressed = -1;// if left -1 and minExpression>0, this will become all samples
	static int minExpression = -1; //if left -1 does nothing
	
	public static void main(String[] args) throws IOException, NotConvergedException, InterruptedException
	{
		String expFile = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/" + "TESTexpression.txt";
//		String expFile = "E:/Groningen/Data/Monogenetic_disease_samples_RadBoud/CountFiles/" + "CountsGENES_Radboud.txt";
//		String expFile = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/" + "est_counts_nocancernocelllineFirst100.txt";
		String chromLocationsFile = "E:/Groningen/Data/GenePositionInfo_23X_24Y_25MT_26rest.txt";

		String writeFolder = expFile.replace(".txt", "/");
		System.out.println(System.getProperty("user.dir"));
		
		boolean writeAll = true;
		boolean correctInputForSTdevs = false;
		boolean correctInputForSTdevsAfterCenter = false;
		boolean log2 = false;
		boolean skipQuantileNorm = true;	
		boolean STdevCutoff = false;
		boolean zScores = false;
		boolean directPCA = false;
		boolean correlation = true;
		
		double correctTotalReadCount = -1;//log((gene+0.5)/total*value)
		double rLog = -1;
		double topVariance = 1;
		double randomValue = 0;
		double duplicateCutoff = 1;
		double highestExpressed = 1;//1 = all genes, 0.5 = 50% highest expressed genes only (removes 50% lowest expressed genes after quantile normalization (then re-normalizes)).
		double spearman = -1;//if 0, does spearman, if >0, it sets all values below this value to 0.
		int ignoreLowestValues = -1;//I still need to fix this
		
		if(args.length == 0)
			checkArgs(args);
		
		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase()){
				case "filename":
					expFile = value;
					break;
				case "chrom":
					chromLocationsFile = value;
					break;
				case "writeall":
					writeAll = Boolean.parseBoolean(value);
					break;
				case "correctstdevs":
					correctInputForSTdevs = Boolean.parseBoolean(value);
					break;
				case "correctstdevsaftercenter":
					correctInputForSTdevsAfterCenter = Boolean.parseBoolean(value);
					break;
				case "correcttotalreadcount":
					correctTotalReadCount = Double.parseDouble(value);
					break;
				case "rlog":
					rLog = Double.parseDouble(value);
					break;
				case "log2":
					log2 = Boolean.parseBoolean(value);
					break;
				case "randomvalue":
					randomValue = Double.parseDouble(value);
					break;
				case "writefolder":
					writeFolder = value;
					break;
				case "duplicate":
					duplicateCutoff = Double.parseDouble(value);
					break;
				case "highestexpressed":
					highestExpressed = Double.parseDouble(value);
					break;
				case "noqn":
					skipQuantileNorm = Boolean.parseBoolean(value);
					break;
				case "stdevcutoff":
					STdevCutoff = Boolean.parseBoolean(value);
					break;
				case "zscores":
					zScores = Boolean.parseBoolean(value);
					break;
				case "topvariance":
					topVariance = Double.parseDouble(value);
					break;
				case "directpca":
					directPCA = Boolean.parseBoolean(value);
					break;
				case "spearman":
					spearman = Double.parseDouble(value);
					break;
				case "correlation":
					correlation = Boolean.parseBoolean(value);
					break;
				case "lowesttoaverage":
					setLowestToAverage = Boolean.parseBoolean(value);
					break;
				case "adjustsampleaverages":
					adjustSampleAverages = Boolean.parseBoolean(value);
					break;	
				case "removegene":
					removeGene = value;
					break;
				case "minsamplesexpressed":
					minSamplesExpressed = Integer.parseInt(value);
					break;
				case "minexpression":
					minExpression = Integer.parseInt(value);
					break;
				default:
					checkArgs(args);
					System.out.println("Incorrect argument supplied; exiting");
					System.exit(1);
			}
		}

		writeParameters(expFile, writeFolder, chromLocationsFile, writeAll, log2, correctInputForSTdevs, 
				randomValue, duplicateCutoff, highestExpressed,skipQuantileNorm, STdevCutoff, zScores, directPCA, spearman,
				ignoreLowestValues, setLowestToAverage, correctTotalReadCount);
		
		run(expFile, writeFolder, chromLocationsFile, writeAll, log2, correctInputForSTdevs,
				randomValue, duplicateCutoff, highestExpressed, skipQuantileNorm, STdevCutoff, correctInputForSTdevsAfterCenter,
				zScores, correctTotalReadCount, topVariance, directPCA, rLog, spearman, correlation, ignoreLowestValues);
	}
	private static void writeParameters(String expFile, String writeFolder, String chromLocationsFile, boolean writeAll,
			boolean log2, boolean correctInputForSTdevs, double randomValue, double duplicateCutoff, double highestExpressed,
			boolean tpm, boolean STdevCutoff, boolean zScores, boolean directPCA, double spearman, int ignoreLowestValues,
			boolean setLowestToAverage, double correctTotalReadCount) throws IOException {
		makeFolder(writeFolder);
		DateFormat dateFormat = new SimpleDateFormat("dd/MM/yyyy HH:mm:ss");
		Date date = new Date();
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(writeFolder+"parameters.txt")));
		writer.write("Date\t" + dateFormat.format(date) + "\n");
		writer.write("Input file\t" + expFile + "\n");
		writer.write("writeFolder\t" + writeFolder + "\n");
		writer.write("chromLocationsFile\t" + chromLocationsFile + "\n");
		writer.write("writeAll\t" + writeAll + "\n");
		writer.write("correctTotalReadCount\t" + correctTotalReadCount + "\n");
		writer.write("log2\t" + log2 + "\n");
		writer.write("correctInputForSTdevs\t" + correctInputForSTdevs + "\n");
		writer.write("randomValue\t" + randomValue + "\n");
		writer.write("duplicateCutoff\t" + duplicateCutoff + "\n");
		writer.write("TopXhighestExpressed\t" + highestExpressed + "\n");
		writer.write("noQN\t" + tpm + "\n");
		writer.write("STdevCutoff\t" + STdevCutoff + "\n");
		writer.write("zScores\t" + zScores + "\n");
		writer.write("directPCA\t" + directPCA + "\n");
		writer.write("spearman\t" + spearman + "\n");
		writer.write("setLowestToAverage\t" + setLowestToAverage + "\n");
		writer.write("adjustSampleAverages\t" + adjustSampleAverages + "\n");
		writer.write("removeGene\t" + removeGene + "\n");
		writer.write("minExpression\t" + minExpression + "\n");
		writer.write("minSamplesExpressed\t" + minSamplesExpressed + "\n");
		writer.close();
	}
	static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
			return;
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
	public static void run(String expFile, String writeFolder, String chromLocationsFile, boolean writeAll, 
			boolean log2, boolean correctInputForSTdevs, double randomAddValue, double removeDuplicates, 
			double highestExpressed, boolean skipQuantileNorm, boolean STdevCutoff, boolean correctInputForSTdevsAfterCenter,
			boolean zScores, double correctTotalReadCount, double topVariance, boolean directPCA, double rLog,
			double spearman, boolean correlation, int ignoreLowestValues) throws IOException, NotConvergedException, InterruptedException
	{		
		pca.PCA.log(" 1. Reading expression file");
		MatrixStruct expressionStruct = new MatrixStruct(expFile);
		
		if(writeAll)
			expressionStruct.write(writeFolder+"startMatrix.txt");
		
		RemoveDuplicates.removeDuplicates(expressionStruct, removeDuplicates, writeFolder,writeAll);

		pca.PCA.log(" 3. Transposing");
		expressionStruct.putGenesOnRows();
		
		//sorting on chromosome locations
		expressionStruct = SortChromosome.sort(expressionStruct, chromLocationsFile);
		
		//adding random values
		expressionStruct.addRandomValues(randomAddValue);
		expressionStruct.write(writeFolder+"randAddedToBelow_" +randomAddValue + ".txt");
		
		if(minExpression > 0)
			expressionStruct.removeLowExpression((writeFolder+"RemovedBelow" +minExpression+ "inAtLeast" + minSamplesExpressed + "samples.txt"), minSamplesExpressed, minExpression);
		
		if(removeGene != null)
			expressionStruct.removeRow(removeGene);
		
		if(!skipQuantileNorm)//quantile Normalize
			QuantileNormalize.quantileNormalize(expressionStruct, writeFolder, writeAll);
		
		if(correctTotalReadCount > 0)
			CorrectReadcounts.correct(writeFolder, correctTotalReadCount, expressionStruct, writeAll, 0.5);
		if(rLog > 0)//I still need to test this function
			RLog.rLog(writeFolder, expressionStruct, rLog, writeAll);
		
		double addVal = 0;
		if(correctTotalReadCount <= 0 && rLog <= 0) // need to add 1 before log to avoid log(0)
			addVal = 1;
		if(log2)
			LogTransform.log2(writeFolder, expressionStruct, writeAll,addVal);
		
		
		
		if(highestExpressed >0 && highestExpressed < 1)
			HighestExpressed.highestExpressed(expressionStruct, skipQuantileNorm, correctTotalReadCount, writeFolder, highestExpressed, STdevCutoff, writeAll);
		
		pca.PCA.log("11. Transposing");
		expressionStruct.transpose();
		
		boolean correl  = true;
		if(topVariance > 0 && topVariance !=1)
			KeepTopVariance.keepTopVariance(expressionStruct, correl, ignoreLowestValues, writeFolder, topVariance);
		
		pca.PCA.log("12 Calculating STdevs");
		System.gc();System.gc();
		MatrixStruct stDevs = expressionStruct.stDevCols();
		stDevs.write(writeFolder + "gene_STDevs.txt");
		System.out.println("Rows = " + expressionStruct.rows()+ " cols = " + expressionStruct.cols());
		
		if(correctInputForSTdevs)
		{
			pca.PCA.log("13 Divide all gene values by STdev for each gene ");	
			expressionStruct.divideBy(stDevs,false);//false corrects columns, true corrects rows
			
			pca.PCA.log("14 Writing matrix divided by gene STdevs");
			expressionStruct.write(writeFolder + "_DividedBySGenesSTdev.txt");
		}	
		
		//expressionStruct.isExpressed(null, expressionStruct.cols()*0.05,10);//genes need to have an expression of at least 1 or more in at least 5% of all the samples.
		if(spearman >= 0)
		{
			Spearman.ranks(writeFolder, expressionStruct, spearman);
			correl = true;
		}
		
		if(setLowestToAverage)//this will cause the lowest values not to contribute to correlation or covariance
			LowestToAverage.lowestToAverage(expressionStruct);
		
		pca.PCA.log("15. Calculating column averages");
		MatrixStruct colAverages = expressionStruct.getAveragesPerCol();
		colAverages.write(writeFolder+ "SAMPLE_Norm_columnAverages.txt");
		pca.PCA.log("16. Centering: Adjusting for column averages");
		expressionStruct.adjustForAverageAllCols(colAverages);
		expressionStruct.write(writeFolder +"SAMPLE_adjustedForGeneAverages.txt");
		pca.PCA.log("17. Calculating row averages");
		MatrixStruct rowAverages = expressionStruct.getAveragesPerRow();
		rowAverages.write(writeFolder+ "SAMPLE_Norm_rowAverages.txt");
		if(adjustSampleAverages)
		{
			pca.PCA.log("18. Adjusting for row averages");
			expressionStruct.adjustForAverageAllrows(rowAverages);
		}
	
		String expNormLogCentFile = writeFolder+"MATRIX_Centered.txt";
		pca.PCA.log("19. Writing centered file in: " + expNormLogCentFile);
		expressionStruct.write(expNormLogCentFile);
		
		double[] cutoffs = KeepTopVariance.getCutoffs(ignoreLowestValues, expressionStruct);
		MatrixStruct geneEigenVectors = null;
		
		if(correctInputForSTdevsAfterCenter)
		{
			MatrixStruct stDevsRows = expressionStruct.stDevRows();
			stDevsRows.write(writeFolder + "_SampleStDevs.txt");
			pca.PCA.log("13 Divide all gene values by STdev for each sample");	
			expressionStruct.divideBy(stDevsRows,true);//false corrects columns, true corrects rows
			
			pca.PCA.log("14 Writing matrix divided by gene STdevs");
			expressionStruct.write(writeFolder + "_DividedBySGenesSTdev.txt");
		}
		
		if(directPCA)
			geneEigenVectors = directPCA(expressionStruct, writeFolder, correlation, null);
		else
			geneEigenVectors = inDirectPCA(expressionStruct, writeFolder, correlation, null);

		pca.PCA.log("25. Calculating PCscores for all samples");
		MatrixStruct[] scoreResults = PCA.scores(geneEigenVectors,expressionStruct, writeFolder+"SAMPLE_PC.txt");
		MatrixStruct PCsampleScores = scoreResults[0];
		
		if(zScores == true)
		{
			pca.PCA.log("26. Calculating Z-scores for all PCscores for all samples");
			MatrixStruct zScoreStats = Zscore.changeToZscores(PCsampleScores);
			
			pca.PCA.log("27. Writing Z-scores");
			PCsampleScores.write(writeFolder+ "pcZscoresSamples.txt");
			zScoreStats.write(writeFolder+ "pcZscores_Stats.txt");
		}
		pca.PCA.log("Files written to: " + writeFolder);
	}
	
	private static MatrixStruct inDirectPCA(MatrixStruct expressionStruct,
			String writeFolder, boolean correlation, double[] cutoffs) throws IOException, NotConvergedException {
		if(expressionStruct.getRowHeaders()[0].contains("ENSG0") && expressionStruct.getRowHeaders()[0].contains("ENST0"))
		{
			System.out.println("Genes should be on Cols and are not, transposing (assumes geneIDs start with ENSG0 or ENST0");
			expressionStruct.transpose();
		}
		String type = "covariance";
		if(correlation)
		{
			expressionStruct.removeNoVariance(writeFolder+"noVarRemoved.txt");
			type = "correlation";
		}
		pca.PCA.log("20. creating covariance matrix over the samples");
		ConcurrentCovariation calculator = new ConcurrentCovariation(20);
		double[][] inMat = expressionStruct.getMatrix();
		double[][] covMatrix = calculator.pairwiseCovariation(inMat,false, null, expressionStruct.getRowHeaders(),correlation, cutoffs);
		
		MatrixStruct covMat = new MatrixStruct(expressionStruct.getRowHeaders(), expressionStruct.getRowHeaders(), covMatrix);
		
		
		pca.PCA.log("21. Writing covariance matrix over the samples");
		String covMatFN = writeFolder+"SAMPLE_"+type+".txt";
		covMat.write(covMatFN);
		
		//to test something
		Matrix covMatOtherFormat = new Matrix(expressionStruct.getRowHeaders(), expressionStruct.getRowHeaders(), covMatrix);
		covMatOtherFormat.getAverageCols(true)
		  .write(covMatFN.replace(".txt", "_Absolute_averages.txt"));
		covMatOtherFormat.write(covMatFN);
		
		inMat = null; covMatrix = null; System.gc();System.gc();
				
		pca.PCA.log("22. calculating eigenvalues over the samples");
		MatrixStruct[] evds = pca.PCA.evd(covMat, Paths.get(writeFolder+"SAMPLE"));
		MatrixStruct eigenVectors = evds[0];
		MatrixStruct PCeigenvalues = evds[1];
		
		pca.PCA.log("23. calculating PCscores over the genes");
		String saveNamePCscoresGene = writeFolder + "GENE_PC.txt";
		MatrixStruct[] PCscoresGenesAndAverages = PCA.scores(eigenVectors,expressionStruct, saveNamePCscoresGene);
		MatrixStruct PCscoresGenes = PCscoresGenesAndAverages[0];
		PCscoresGenesAndAverages=null;
		System.gc();
		
		pca.PCA.log("24. Transform PCscores to eigenvectors of Genes");
		String saveNameEigenVectorsOverGenes = writeFolder + "GENE.eigenvectors.txt";
		MatrixStruct geneEigenVectors = PCA.transform(PCscoresGenes, PCeigenvalues, saveNameEigenVectorsOverGenes);
		return geneEigenVectors;
	}
	private static MatrixStruct directPCA(MatrixStruct expressionStruct, String writeFolder, boolean correlation, double[] cutoffs) throws IOException, InterruptedException {
		if(!expressionStruct.getRowHeaders()[0].contains("ENSG0") && !expressionStruct.getRowHeaders()[0].contains("ENST0"))
		{
			System.out.println("Genes should be on rows and are not, transposing (assumes geneIDs start with ENSG0 or ENST0");
			expressionStruct.transpose();
		}
		double[][] inMat = expressionStruct.getMatrix();
		ConcurrentCovariation calculatorGenes = new ConcurrentCovariation(20);
		String type = "covariance";
		if(correlation)
		{
			expressionStruct.removeNoVariance(writeFolder+"noVarRemoved.txt");
			type = "correlation";
		}
		String covMatFN = writeFolder+"gene_"+type+".txt";
		double[][] covMatrix = calculatorGenes.pairwiseCovariation(inMat,false, null, expressionStruct.getRowHeaders(),correlation, cutoffs);//last argument, if false = covariance; true = correlation.
		
		Matrix covMat = new Matrix(expressionStruct.getRowHeaders(), expressionStruct.getRowHeaders(), covMatrix);
		pca.PCA.log("21. Writing covariance matrix over the samples");
		//get and write averages per Gene
		covMat.getAverageCols(true)
			  .write(covMatFN.replace(".txt", "_Absoulte_averages.txt"));
		
		covMat.write(covMatFN);
		inMat = null; covMatrix = null; System.gc();System.gc();
		
		inMat = null; System.gc();System.gc();
		
		PCA.log("Matrix decomposition");
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File("export.sh")));
		//export DYLD_LIBRARY_PATH="/opt/intel/mkl/lib/":"/opt/intel/compilers_and_libraries_2016.1.111/mac/compiler/lib":$DYLD_LIBRARY_PATH
		writer.write("export DYLD_LIBRARY_PATH=\"/opt/intel/mkl/lib/\":\"/opt/intel/compilers_and_libraries_2016.1.111/mac/compiler/lib\":$DYLD_LIBRARY_PATH\n");
		writer.write("/Volumes/Promise_RAID/juha/PCA++/pca evd " + covMatFN +"\n");
		writer.write("mv eigenvectors.txt " + writeFolder + "GENE.eigenvectors.txt"+"\n");
		writer.write("mv eigenvalues.txt " + writeFolder + "GENE.eigenvalues.txt"+"\n");
		writer.close();		
		String command = "sh export.sh";
		Process p;
		p = Runtime.getRuntime().exec(command);
		p.waitFor();
		MatrixStruct geneEigenVectors = new MatrixStruct(writeFolder + "GENE.eigenvectors.txt",-1,11000);//only use the first 11000 eigenvectors
		expressionStruct.transpose();
		return geneEigenVectors;
	}
	
	static void makeFolder(String writeFolder) 
	{
		File folder = new File(writeFolder);
		if(!folder.exists())
		{
			folder.mkdir();
		}
		
	}
}
