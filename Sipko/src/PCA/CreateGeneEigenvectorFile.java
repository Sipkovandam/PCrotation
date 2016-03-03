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
import java.util.Date;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

import no.uib.cipr.matrix.NotConvergedException;
import pca.MatrixStruct;
import pca.PCA;

public class CreateGeneEigenvectorFile 
{
	public static void main(String[] args) throws IOException, NotConvergedException
	{
		String expFile = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/" + "TESTexpression.txt";
		//String expFile = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/" + "4DownSyndrome3Normal3Cancer_counts.txt";
		String chromLocationsFile = "E:/Groningen/Data/GenePositionInfo_23X_24Y_25MT_26rest.txt";

		String writeFolder = expFile.replace(".txt", "/");
		System.out.println(System.getProperty("user.dir"));
		
		boolean writeAll = true;
		boolean correctInputForSTdevs = false;
		boolean log2 = true;
		boolean tpm = false;	
		boolean STdevCutoff = false;
		
		double randomValue = 0;
		double duplicateCutoff = 1;
		double highestExpressed = 1;//1 = all genes, 0.5 = 50% highest expressed genes only (removes 50% lowest expressed genes after quantile normalization (then re-normalizes)).
		
		
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
				case "tpm":
					tpm = Boolean.parseBoolean(value);
					break;
				case "stdevcutoff":
					STdevCutoff = Boolean.parseBoolean(value);
					break;
				default:
					checkArgs(args);
					System.out.println("Incorrect argument supplied; exiting");
					System.exit(1);
			}
		}

		writeParameters(expFile, writeFolder, chromLocationsFile, writeAll, log2, correctInputForSTdevs, randomValue, duplicateCutoff, highestExpressed,tpm, STdevCutoff);
		
		run(expFile, writeFolder, chromLocationsFile, writeAll, log2, correctInputForSTdevs, randomValue, duplicateCutoff, highestExpressed, tpm, STdevCutoff);
	}
	private static void writeParameters(String expFile, String writeFolder, String chromLocationsFile, boolean writeAll,
			boolean log2, boolean correctInputForSTdevs, double randomValue, double duplicateCutoff, double highestExpressed,
			boolean tpm, boolean STdevCutoff) throws IOException {
		makeFolder(writeFolder);
		DateFormat dateFormat = new SimpleDateFormat("dd/MM/yyyy HH:mm:ss");
		Date date = new Date();
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(writeFolder+"parameters.txt")));
		writer.write("Date\t" + dateFormat.format(date) + "\n");
		writer.write("Input file\t" + expFile + "\n");
		writer.write("writeFolder\t" + writeFolder + "\n");
		writer.write("chromLocationsFile\t" + chromLocationsFile + "\n");
		writer.write("writeAll\t" + writeAll + "\n");
		writer.write("log2\t" + log2 + "\n");
		writer.write("correctInputForSTdevs\t" + correctInputForSTdevs + "\n");
		writer.write("randomValue\t" + randomValue + "\n");
		writer.write("duplicateCutoff\t" + duplicateCutoff + "\n");
		writer.write("TopXhighestExpressed\t" + highestExpressed + "\n");
		writer.write("tpm\t" + tpm + "\n");
		writer.write("STdevCutoff\t" + STdevCutoff + "\n");
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
				+ "tpm=<false/true> - whether expression values is tpm normalized (if so no Quantile normalization is applied) (default=false)\n"
				+ "highestExpressed=<number> - Percentage of genes to keep (genes most highly expressed on average)");
		System.exit(1);
	}
	public static void run(String expFile, String writeFolder, String chromLocationsFile, boolean writeAll, 
			boolean log2, boolean correctInputForSTdevs, double randomAddValue, double removeDuplicates, double highestExpressed, boolean tpm, boolean STdevCutoff) throws IOException, NotConvergedException
	{		
		pca.PCA.log(" 1. Reading expression file");
		MatrixStruct expressionMatrixStruct = new MatrixStruct(expFile);
		
		pca.PCA.log(" 2. Removing duplicates (r>"+removeDuplicates+")");
//		expressionMatrixStruct.write(writeFolder+"startMatrix.txt");
		expressionMatrixStruct = removeDuplicates(expressionMatrixStruct, removeDuplicates);
		String duplicatesRemovedFN = writeFolder+"DuplicatesRemoved.txt";
		if(writeAll)expressionMatrixStruct.write(duplicatesRemovedFN);
		
		pca.PCA.log(" 3. Transposing");
		expressionMatrixStruct.transpose();
			
		if(chromLocationsFile != null)
		{
			pca.PCA.log("Sorting IDs based on chromosome locations");
			MatrixStruct chromLocs = new MatrixStruct(chromLocationsFile);
			expressionMatrixStruct = SortChromosome.sort(expressionMatrixStruct, chromLocs);
		}
	
		if(randomAddValue > 0)
		{
			pca.PCA.log(" 4. Adding random values <"+ randomAddValue +" to values <" + randomAddValue);
			expressionMatrixStruct.addRandomValues(randomAddValue);
			expressionMatrixStruct.write(writeFolder+"randAddedToBelow_" +randomAddValue + ".txt");
		}

		Matrix expressionMatrix = new Matrix(expressionMatrixStruct);
		expressionMatrixStruct = null;
			
		pca.PCA.log(" 5. Removing genes without variance");
		expressionMatrix.removeNoVariance();
		
		pca.PCA.log(" 5.1 Writing file from which genes without variance are removed");
		String removedGenesFN = writeFolder+"noVarRemoved.txt";
		if(writeAll)expressionMatrix.write(removedGenesFN);
		
		if(highestExpressed >0 && highestExpressed < 1)
		{
			if(!tpm)
			{
				pca.PCA.log(" Quantile normalization before taking averageCutoff");
				expressionMatrix.quantileNormAdjust(expressionMatrix.quantileNormVector());
			}	
			String averagesFN = writeFolder+"NormalizedAveragesAllGenes.txt";
			expressionMatrix.calcAvgRows()
							.write(averagesFN);
			
			expressionMatrixStruct = new MatrixStruct(expressionMatrix.rowNames, expressionMatrix.colNames, expressionMatrix.values);
			pca.PCA.log("  . Removing the " + ((1.0-highestExpressed)*100) + " percent lowest expressed genes" );
			
			if(STdevCutoff)
			{
				MatrixStruct stDevs = expressionMatrixStruct.stDevRows();
				String stdevFN = writeFolder + "gene_STDevsForCutoff.txt";
				stDevs.write(stdevFN);
				PCcorrection.keepTopPercentage(expressionMatrixStruct,stdevFN, highestExpressed, averagesFN.replace(".txt", "top_" + highestExpressed+ ".txt"), true);
			}
			else
				PCcorrection.keepTopPercentage(expressionMatrixStruct,averagesFN, highestExpressed, averagesFN.replace(".txt", "top_" + highestExpressed+ ".txt"), false);
			
			expressionMatrix = new Matrix(expressionMatrixStruct);
		}
		
		if(!tpm)
		{
			pca.PCA.log(" 6. Calculating quantile normalization vector");
			Matrix qNormVector = expressionMatrix.quantileNormVector();
			if(writeAll)qNormVector.write(writeFolder+ "SAMPLE_QuantileVector.txt");
		
		
			pca.PCA.log(" 7. Quantile normalization");
			expressionMatrix.quantileNormAdjust(qNormVector);
			String quantFNnotLogged = writeFolder+ "SAMPLE_QuantileNormalized.txt";
			pca.PCA.log(" 8. Writing quantile normalized data in: " + quantFNnotLogged);
			if(writeAll)expressionMatrix.write(quantFNnotLogged);
		}
		
		if(log2)
		{
			pca.PCA.log(" 9. Log2 transforming");
			expressionMatrix.log2Transform();
			String quantFN = writeFolder+ "SAMPLE_NormalizedLog2.txt";
			
			pca.PCA.log("10. Writing logged SAMPLE_Log2 normalized data in: " + quantFN);
			if(writeAll)expressionMatrix.write(quantFN);
		}
		
		pca.PCA.log("11. Transposing");
		expressionMatrix.transpose();
		
		pca.PCA.log("12 Calculating STdevs");		
		MatrixStruct expressionStruct = new MatrixStruct(expressionMatrix.rowNames, expressionMatrix.colNames, expressionMatrix.values); //I still need to rewrite add a few functions in Matrix to MatrixStruct, but cba atm
		System.gc();System.gc();
		MatrixStruct stDevs = expressionStruct.stDevCols();
		stDevs.write(writeFolder + "gene_STDevs.txt");
		if(correctInputForSTdevs)
		{
			pca.PCA.log("13 Divide all gene values by STdev for each gene");	
			expressionStruct.divideBy(stDevs,false);//false corrects columns, true corrects rows
			
			pca.PCA.log("14 Writing matrix divided by gene STdevs");
			expressionStruct.write(writeFolder + "_DividedBySGenesSTdev.txt");
			expressionMatrix = new Matrix(expressionStruct);
		}	

		pca.PCA.log("15. Calculating column averages");
		Matrix colAverages = expressionMatrix.calcAvgCols();//NOTE: this is a vector that is based on rows already being corrected for averages.
		colAverages.write(writeFolder+ "SAMPLE_Norm_columnAverages.txt");
		
		pca.PCA.log("16. Centering: Adjusting for column averages");
		expressionMatrix.adjustForAverageAllCols(colAverages);
		expressionMatrix.write(writeFolder +"SAMPLE_adjustedForGeneAverages.txt");
		
		pca.PCA.log("17. Calculating row averages");
		Matrix rowAverages = expressionMatrix.calcAvgRows();
		rowAverages.write(writeFolder+ "SAMPLE_Norm_rowAverages.txt");
		
		pca.PCA.log("18. Adjusting for row averages");
		expressionMatrix.adjustForAverageAllrows(rowAverages);
		
		String expNormLogCentFile = writeFolder+"MATRIX_Centered.txt";
		pca.PCA.log("19. Writing centered file in: " + expNormLogCentFile);
		expressionMatrix.write(expNormLogCentFile);
			
		expressionStruct = new MatrixStruct(expressionMatrix.rowNames, expressionMatrix.colNames, expressionMatrix.values);
		expressionMatrix = null;
		
		pca.PCA.log("20. creating covariance matrix over the samples");
		ConcurrentCovariation calculator = new ConcurrentCovariation(20);
		double[][] inMat = expressionStruct.getSquareMatrix();
		double[][] covMatrix = calculator.pairwiseCovariation(inMat);
		MatrixStruct covMat = new MatrixStruct(expressionStruct.getRowHeaders(), expressionStruct.getRowHeaders(), covMatrix);
		
		pca.PCA.log("21. Writing covariance matrix over the samples");
		covMat.write(writeFolder+"SAMPLE_covariance.txt");
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
		
		pca.PCA.log("25. Calculating PCscores for all samples");
		MatrixStruct[] scoreResults = PCA.scores(geneEigenVectors,expressionStruct, writeFolder+"SAMPLE_PC.txt");
		MatrixStruct PCsampleScores = scoreResults[0];
		
		pca.PCA.log("26. Calculating Z-scores for all PCscores for all samples");
		MatrixStruct zScoreStats = Zscore.changeToZscores(PCsampleScores);
		
		pca.PCA.log("27. Writing Z-scores");
		PCsampleScores.write(writeFolder+ "pcZscoresSamples.txt");
		zScoreStats.write(writeFolder+ "pcZscores_Stats.txt");
		
		pca.PCA.log("Files written to: " + writeFolder);
	}
	private static MatrixStruct removeDuplicates(MatrixStruct expression, double duplicateCutoff) {
		//this function assumes duplicate rows are always next to each other (as is usually the case)
		//This saves some computational time
		ArrayList<Integer> rowsToRemove = new ArrayList<Integer>();
		for(int r = 0; r < expression.rows()-1; r++)
		{
			double correlation = Correlation.correlate(expression.getRowValues(r), expression.getRowValues(r+1));
			if(correlation > duplicateCutoff)
				rowsToRemove.add(r+1);
		}
		System.out.println("Removing  " + rowsToRemove.size() + " duplicates");
		MatrixStruct adjustedMatrix = new MatrixStruct(expression.rows()-rowsToRemove.size(),expression.cols());
		adjustedMatrix.setColHeaders(expression.getColHeaders());
		int a = 0;
		int outRow = 0;
		if(rowsToRemove.size()==0)
			return expression;
		for(int r = 0; r < expression.rows(); r++)
		{
			int skip = rowsToRemove.get(a);
			if(skip == r)
			{
				if(a<rowsToRemove.size()-1)
					a++;
				continue;
			}
			adjustedMatrix.setRow(outRow, expression.getRowHeaders()[r], expression.getRowValues(r));
			outRow++;
		}
		
		return adjustedMatrix;
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
