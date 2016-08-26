package PrepData;

import java.io.IOException;
import java.util.Hashtable;
import java.util.List;

import PCA.MatrixStruct;
import PCA.RLog;
import PCA.SortChromosome;

public class PCA_createTestMatrix 
{
	//creates 2 test matrixes, 1 that represents the "public" data and 1 that represents the diseases.
	
	static String publicHealthyFN = "/Volumes/Promise_RAID/GeneNetwork/Sipko/CountFiles/est_counts_nocancernocellline.txt.gz";
	static String clinicFN = "/Volumes/Promise_RAID/sipko/Clinic/Counts/Kallistocounts_GENES.txt.gz";
	static String downFN = "/Volumes/Promise_RAID/GeneNetwork/Sipko/CountFiles/18DownSyndrome26Normal2Cancer_counts.txt.gz";
	static String chromLocationsFN = "/Volumes/Promise_RAID/GeneNetwork/Sipko/GenePositionInfo_23X_24Y_25MT.txt";
	static String writeFolder = "/Volumes/Promise_RAID/GeneNetwork/Sipko/temp/";
	
//	static String publicHealthyFN = "E:/Groningen/Test/PCA_createTestMatrix/TESTexpression.txt";
//	static String clinicFN = "E:/Groningen/Data/RNAseq_clinic/Analyse/CountsAll_22-08-2016/Kallistocounts_GENES.txt.gz";
//	static String downFN = "E:/Groningen/Data/PublicSamples/DownSyndrome/18DownSyndrome26Normal2Cancer_counts.txt";
//	static String chromLocationsFN = "E:/Groningen/Data/GenePositionInfo_23X_24Y_25MT.txt";
//	static String writeFolder = "E:/Groningen/Test/PCA_createTestMatrix/";
	
	static int nGenes = 10000;
	static int nSamples= 2000;
	
	public static void main(String[] args) throws IOException
	{
		if(args.length ==1)
			nGenes=Integer.parseInt(args[0]);
		
		if(args.length ==2)
		{
			nGenes=Integer.parseInt(args[0]);
			nSamples= Integer.parseInt(args[1]);
		}
		MatrixStruct down = new MatrixStruct(downFN);
		down.putGenesOnRows();
		System.out.println("rows1 = " + down.rows());
		
		String writeDiseaseFN = writeFolder+"disease.txt";
		String writeHealthyFN  = writeFolder+"healthy.txt";
		
		MatrixStruct geneAvg= getAverageExpression(down);//also normalizes
		MatrixStruct expressionStruct=keepTopExpressedIDs(geneAvg);
		expressionStruct = getRandomSamples(expressionStruct);
		
		MatrixStruct clinic = new MatrixStruct(clinicFN);
		clinic.putGenesOnRows();
		
		MatrixStruct merge = clinic.mergeColumns(down);//need to fix this
		merge.keepRows(expressionStruct);
		merge = SortChromosome.sort(merge, chromLocationsFN);//order the genes by their chromosome number and location
		expressionStruct = SortChromosome.sort(expressionStruct, chromLocationsFN);//order the genes by their chromosome number and location
		System.out.println(merge.getRowHeaders()[0]);
		System.out.println(expressionStruct.getRowHeaders()[0]);
		merge.write(writeDiseaseFN);
		expressionStruct.write(writeHealthyFN);
	}

	private static MatrixStruct getRandomSamples(MatrixStruct expressionStruct) 
	{
		MatrixStruct subset = new MatrixStruct(expressionStruct.rows(), nSamples);
		subset.setRowHeaders(expressionStruct.getRowHeaders());
		int stepsize = expressionStruct.cols()/nSamples;//to ensure you get samples from as wide a range of backgrounds as possible (so it has the same distribution as the actual data)
		if(stepsize ==0)
			stepsize =1;
		for(int n = 0; n < nSamples; n++)
		{
			int index = stepsize*n;
			//System.out.println(index);
			double[] values = expressionStruct.getColValues(index);
			subset.setColValues(values,n);
			subset.setColHeader(n, expressionStruct.getColHeaders()[index]);
		}
			
		return subset;
	}

	private static MatrixStruct keepTopExpressedIDs(MatrixStruct geneAvg) {
		MatrixStruct expressionStruct = new MatrixStruct(publicHealthyFN);//reload matrix
		expressionStruct.putGenesOnRows();//put genes on rows if they are not already
		
		Hashtable<String,Integer> toKeep = new Hashtable<String,Integer>();
		for(int g = 0; g < nGenes; g++)
		{
			toKeep.put(geneAvg.getRowHeaders()[g], g);
		}
		expressionStruct.keepIDs(toKeep);
		return expressionStruct;
	}

	private static MatrixStruct getAverageExpression(MatrixStruct down) throws IOException {
		MatrixStruct expressionStruct = new MatrixStruct(publicHealthyFN);//load matrix
		expressionStruct.putGenesOnRows();//put genes on rows if they are not already
		expressionStruct.keepRows(down);
		MatrixStruct chrom = new MatrixStruct(chromLocationsFN);
		expressionStruct.keepRows(chrom);
		RLog.rLog(writeFolder, expressionStruct, false, writeFolder+"geomean.txt");
		expressionStruct.logTransform(2, 0.5);
		MatrixStruct geneAvg = expressionStruct.getAveragesPerRow();//get the average expression per gene
		geneAvg.sortCol(0);
		geneAvg.write(writeFolder+"geneAverages.txt");
		//keep only the 10.000 highest expressed genes;
		return geneAvg;
	}
}
