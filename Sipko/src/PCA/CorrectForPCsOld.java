package PCA;

import java.io.IOException;
import java.util.ArrayList;

import pca.MatrixStruct;
import pca.PCA;

public class CorrectForPCsOld {
//depricated class
	public static void main(String[] args) throws IOException 
	{
//		String sampleFile = "E:/Groningen/Data/PublicSamples/Test7/" + "4DownSyndrome3Normal3Cancer_countsRND.txt";
//		String vectorFolder = "E:/Groningen/Data/PublicSamples/Test7/est_counts_nocancernocelllineCorrelationMatrixPCAoverGenes/";
		
		String sampleFile = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/" + "RandomSamples.txt";
		String vectorFolder = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/TESTexpression/";

		String PCsToAdjust = "0";//calling this function with 0 will just rotate the matrix back without adjusting for any PCs.
		
		checkArgs(args);
		if(!System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
		{
			sampleFile = args[0];
			vectorFolder = args[1]+"/";
			PCsToAdjust = args[2];
		}
		
		String writeFolder = sampleFile.replace(".txt", "_Adj/");
		MatrixStruct[] PCscoresSample = RotateSample.rotate(sampleFile, vectorFolder, writeFolder, true,true,false);
		pca.PCA.log("9.1 Calculating zScores");
		String zScoreStats = vectorFolder+"pcZscores_Stats.txt";
		MatrixStruct zScoreMatrix = PCscoresSample[0].copy();
		Zscore.changeToZscores(zScoreMatrix, zScoreStats);
		pca.PCA.log("9.2 Writing zScores");
		zScoreMatrix.write(writeFolder+"pcZscoresSamples.txt");

		/**set PCs to 0 (correcting for principal components)**/
		ArrayList<Integer> PCsToCorrect = parsePCs(PCsToAdjust);
		//int[] PCsToCorrect = new int[]{1,2,3,4,5};
		pca.PCA.log("10. Setting PCs to 0");
		PCscoresSample[0].setZeroRows(PCsToCorrect);
		
		//can also adjust by Names of the PC, but need to know the actual PC names in the file
		//String[] PCsToCorrect = new String[]{"PC1","PC2","PC3","PC4","PC5"};
		//PCscoresSingleSample.setZeroRows(PCsToCorrect);
		
		/**Rotate sample back to original expression**/
		String saveNameSingleRotatedBack = writeFolder + "SAMPLE_RotatedBack.txt";
		MatrixStruct geneEigenVectors = new MatrixStruct(vectorFolder+"GENE.eigenvectors.txt");
		String signalPerPCFolder = writeFolder+ "signalPerPC/";
//		pca.PCA.log("11. Calculating expression for each gene for each PC for each sample");
//		calculateGeneExpressionPerPC(PCscoresSample,geneEigenVectors,signalPerPCFolder);
		pca.PCA.log("12. Rotate sample back to original expression");
		MatrixStruct[] rotatedBack = PCA.rotateBack(geneEigenVectors,PCscoresSample[0], saveNameSingleRotatedBack);
		rotatedBack[0].setColHeaders(PCscoresSample[2].getRowHeaders());
		String rotatedAdjusted = writeFolder+"Adjusted_PC"+PCsToAdjust+".txt";;
		rotatedBack[0].write(rotatedAdjusted);
		rotatedBack[0].transpose();
		rotatedBack[0].write(rotatedAdjusted.replace(".txt","_transposed.txt"));
		pca.PCA.log("Corrected file written to: " + rotatedAdjusted);
//		/**Add averages that were initially removed**/
//		pca.PCA.log("12. Adding column averages again");
//		rotatedBack[0].addAveragesCols(new MatrixStruct(vectorFolder+"SAMPLE_QuantNorm_columnAverages.txt"));
//		pca.PCA.log("13. Adding row averages again");
//		rotatedBack[0].addAveragesRows(new MatrixStruct(writeFolder+"rowAverages.txt"));
//		String writeName = writeFolder+"Adjusted_PC"+PCsToAdjust+"_Uncentered.txt";
//		pca.PCA.log("14. Writing adjusted file to:");
//		rotatedBack[0].write(writeName);
//		pca.PCA.log("File written to: " + writeName);
//		pca.PCA.log("15. Transposing adjusted file");
//		rotatedBack[0].transpose();
//		String transposedWriteName = writeName.replace(".txt", "_transposed.txt"); 
//		pca.PCA.log("16. Writing transposed file to" + transposedWriteName);
//		rotatedBack[0].write(transposedWriteName);
	}

	private static void calculateGeneExpressionPerPC(MatrixStruct[] PCscoresSample, MatrixStruct geneEigenVectors, String signalPerPCFolder) throws IOException {
		RotateSample.makeFolder(signalPerPCFolder);
		for(int sample = 0; sample < PCscoresSample[0].cols();sample++)
		{
			pca.PCA.log("Sample: " + PCscoresSample[0].getColHeaders()[sample] + " " +sample +"/" +PCscoresSample[0].cols());
			MatrixStruct samplePCsignals = new MatrixStruct(geneEigenVectors.cols(), geneEigenVectors.rows());
			samplePCsignals.setRowHeaders(geneEigenVectors.getColHeaders());//geneEigenVectors.colHeaders are genes
			samplePCsignals.setColHeaders(geneEigenVectors.getRowHeaders());//geneEigenVectors.rowHeaders are PCs
			for(int pc = 0; pc < geneEigenVectors.rows(); pc++)
			{
				for(int gene = 0; gene < geneEigenVectors.cols(); gene++)
				{
					if(gene == 0 && pc == 0 && sample == 2)
						System.out.println(PCscoresSample[0].matrix.get(pc, sample)+"*"+ geneEigenVectors.matrix.get(pc,gene) + "=" + PCscoresSample[0].matrix.get(pc, sample)*geneEigenVectors.matrix.get(pc,gene));
					//multiply PCscore * eigenvectorCoefficient for that gene for the current PC
					double geneValue =  PCscoresSample[0].matrix.get(pc, sample)*geneEigenVectors.matrix.get(pc,gene);
					samplePCsignals.matrix.set(gene, pc, geneValue);
				}	
			}
			samplePCsignals.write(signalPerPCFolder+PCscoresSample[0].getColHeaders()[sample]+".txt");
		}
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
}
