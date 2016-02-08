package PCA;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import pca.MatrixStruct;

public class CorrectSampleForPCsLudeApproach 
{
	public static void main(String[] args) throws IOException
	{
		String sampleFile = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/" + "RandomSamples.txt";
		String vectorFolder = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/TESTexpression/";
		
//		String sampleFile = "E:/Groningen/Data/PublicSamples/DownSyndrome/" + "18DownSyndrome26Normal274Cancer_counts2.txt";
//		String vectorFolder = "E:/Groningen/Data/PublicSamples/DownSyndrome/est_counts_nocancernocellline/";
//		String vectorFolder = "E:/Groningen/Data/PublicSamples/ServerTest4/est_counts_nocancernocellline/";

		
		
		String PCsToAdjust = "1-50";//calling this function with 0 will just rotate the matrix back without adjusting for any PCs.
		
		checkArgs(args);
		if(!System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko") && args.length !=0)
		{
			sampleFile = args[0];
			vectorFolder = args[1]+"/";
			PCsToAdjust= args[2];
		}
		String writeFolder = sampleFile.replace(".txt", "_AdjLude/");
		
		/**set PCs to 0 (correcting for principal components)**/
		ArrayList<Integer> PCsToCorrect = parsePCs(PCsToAdjust);
		
		/**6. Calculate PCscores for single sample**/
		makeFolder(writeFolder);
		pca.PCA.log(" 1. Loading sample matrix");
		Matrix singleSample = new Matrix(sampleFile);//expressionMatrix.getRow(0);
		
		pca.PCA.log(" 2. Transposing");
		singleSample.transpose();
		Matrix quantVector = new Matrix(vectorFolder+"SAMPLE_QuantileVector.txt");
		pca.PCA.log(" 3. Removing rows that do not exist in the quantile normalization vector");
		singleSample.keepRows(quantVector);
		String rowsRemovedFileName = writeFolder+"_rowsNotInQuantVectorRemoved.txt";
		singleSample.write(rowsRemovedFileName);
		pca.PCA.log(" 4. Quantile normalization adjustion");
		singleSample.quantileNormAdjust(quantVector);
		String quantileAdjustedFN = writeFolder+ "Quantile_adjusted.txt";
		pca.PCA.log(" 5. Writing quantile normalization adjusted file to:" + quantileAdjustedFN);
		pca.PCA.log("Log transforming");
		singleSample.write(quantileAdjustedFN.replace(".txt", "NotLogged.txt"));
		singleSample.log2Transform();//Doing this after the quantile normalization now
		singleSample.write(quantileAdjustedFN);
		singleSample.transpose();
		pca.PCA.log(" 6. Adjusting for column averages (centering to target PC space)");
		singleSample.adjustForAverageAllCols(new Matrix(vectorFolder+"SAMPLE_QuantNorm_columnAverages.txt"));
		String writeColsCentered = writeFolder+"ColsCentered.txt";
		singleSample.write(writeColsCentered);
		//pca.PCA.log(" 7. Adjusting for row averages (centering to target PC space)");
		Matrix averages = singleSample.calcAvgRows();
		String rowAveragesFileName = writeFolder+"rowAverages.txt";
		averages.write(rowAveragesFileName);
		//singleSample.adjustForAverageAllrows(averages);<-- I tested both with and without, it makes no difference (biggest value change was 10-7).
		String centeredFN = writeFolder+ "Quantile_adjusted.centered.txt";
		pca.PCA.log(" 8. Writing PC centered file to: " + centeredFN);
		singleSample.write(centeredFN);		
		pca.PCA.log(" 9. Reading eigenvectors over the genes: ");
		MatrixStruct geneEigenVectors = new MatrixStruct(vectorFolder+"GENE.eigenvectors.txt");
		//MatrixStruct geneEigenVectors = new MatrixStruct(vectorFolder+"GENE.eigenvectorsFirst300.txt");
		pca.PCA.log(" 9. Regressing out signals: ");
		for(int sample = 0; sample< singleSample.rowNames.length;sample++)
		{
			for(int pc : PCsToCorrect)
			{
				double[] rowEigenVector = geneEigenVectors.getRowValues(pc);
			
				double[] rc = getLinearRegressionCoefficients(rowEigenVector,singleSample.values[sample]);
				for (int gene=0; gene<singleSample.colNames.length; gene++) 
				{
					//if(gene < 10 && pc<5)
					//	System.out.println("gene " + singleSample.colNames[gene] + " rc[0]:" + rc[0] + " rc1: "+ rc[1] + " pc: " + pc + " rowEigenVector[gene]*rc[0]:" + rowEigenVector[gene]*rc[0]);
					singleSample.values[sample][gene]-=rowEigenVector[gene]*rc[0] + rc[1];
				}
			}
		}
		String regressedOutFileName = writeFolder+"AdjustedFor" + PCsToAdjust + ".txt";
		pca.PCA.log("10. Writing corrected file to file to: " + regressedOutFileName);
		singleSample.write(regressedOutFileName);
		singleSample.transpose();
		singleSample.write(regressedOutFileName.replace(".txt", "_transposed.txt"));
		
//		MatrixStruct regressedSample = new MatrixStruct(singleSample.rowNames,singleSample.colNames,singleSample.values);
//		/**Add averages that were initially removed**/
//		pca.PCA.log("12. Adding column averages again");
//		regressedSample.addAveragesCols(new MatrixStruct(vectorFolder+"SAMPLE_QuantNorm_columnAverages.txt"));
//		pca.PCA.log("13. Adding row averages again");
//		regressedSample.addAveragesRows(new MatrixStruct(rowAveragesFileName));
//		String writeName = regressedOutFileName.replace(".txt","_averagesReAdded.txt");
//		pca.PCA.log("14. Writing adjusted file to:");
//		regressedSample.write(writeName);
//		pca.PCA.log("File written to: " + writeName);
//		pca.PCA.log("15. Transposing adjusted file");
//		regressedSample.transpose();
//		String transposedWriteName = writeName.replace(".txt", "_transposed.txt"); 
//		pca.PCA.log("16. Writing transposed file to" + transposedWriteName);
//		regressedSample.write(transposedWriteName);
	}
	public static ArrayList<Integer> parsePCs(String PCsToAdjust) 
	{//function takes format like: "1,4,5-10,3-10"
		ArrayList<Integer> PCs = new ArrayList<Integer>();
		if(PCsToAdjust.contains("-"))
		{
			String[] eles = PCsToAdjust.split("-");
			for(int e = 0; e < eles.length-1; e++)
			{
				String[] ele = eles[e].split(",");
				int last = ele.length-1;
				int start = Integer.parseInt(ele[last]);
				String[] ele2 = eles[e+1].split(",");
				int end = Integer.parseInt(ele2[0]);
				for(int n = start; n <= end; n++)
					PCs.add(n-1);
			}
		}
		if(PCsToAdjust.contains(","))
		{
			String[] eles = PCsToAdjust.split(",");
			for(int e = 0; e < eles.length; e++)
			{
				if(!eles[e].contains("-"))
					PCs.add(Integer.parseInt(eles[e])-1);
			}
		}
		
//		System.out.println("size = " + PCs.size());
//		for(int PC : PCs)
//		{
//			System.out.println("pc = " + PC);
//		}
		
		return PCs;
	}
	public static void checkArgs(String[] args)
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
			return;
		if(args.length != 3)
		{
			System.out.println("This script removes the signal from 1 or multiple PCs from a sample.\n"
					+ "Note that the script will not return to orginial values, but to quantile normalized values\n"
					+ "It uses the following 3 arguments:\n"
					+ "1. expression File of the samples \n"
					+ "2. Folder wehre the vector and corresponding files, required for the rotation ar elocated\n"
					+ "   This folder can be created using CreateGeneEigenvectorFile.java script\n"
					+ "3. PCs to correct for. The following formats can be used:"
					+ "   1-100"
					+ "   1,2,6,8"
					+ "   or a combination: 1,5-10,66,100-200");
			System.exit(1);
		}
	}
	private static void makeFolder(String writeFolder) 
	{
		File folder = new File(writeFolder);
		if(!folder.exists())
		{
			folder.mkdir();
		}
		
	}

	public static double[] getLinearRegressionCoefficients(double[] xVal, double[] yVal) {
        double n = (double) xVal.length;
        double sumX = 0; double sumXX = 0;
        double sumY = 0;
        double sumXY = 0;
        for (int x=0; x<xVal.length; x++) {
            sumX+= xVal[x];
            sumXX+=xVal[x] * xVal[x];
            sumY+= yVal[x];
            sumXY+=xVal[x] * yVal[x];
        }
        double sXX = sumXX - sumX*sumX / n;
        double sXY = sumXY - sumX*sumY / n;
        double a = sXY / sXX;
        double b = (sumY - a * sumX) / n;
        double[] regressionCoefficients = new double[2];
        if(Double.isNaN(a) || Double.isNaN(b))
        {
        	a = 0;
        	b = 0;
        }
        regressionCoefficients[0] = a;
        regressionCoefficients[1] = b;
        return regressionCoefficients;
    }
}
