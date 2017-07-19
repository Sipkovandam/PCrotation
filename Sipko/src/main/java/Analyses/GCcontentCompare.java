package Analyses;

import java.io.IOException;
import java.util.Enumeration;
import java.util.Hashtable;

import MatrixScripts.MatrixStruct;

public class GCcontentCompare 
{
	public static String trueGC_FN = "E:/Groningen/Data/GCcontent/GCcontentPublicSamples.txt";
	public static String estimatedGC_FN = "E:/Groningen/Data/GCcontent/est_counts_22214_samples_transposedGENES_GCcontent.txt";
	public static String writeFN = "E:/Groningen/Data/GCcontent/GCcompare_result.txt";
	
	public static void main (String[] args) throws IOException
	{
		MatrixStruct trueGC = new MatrixStruct(trueGC_FN);
		MatrixStruct estimatedGC = new MatrixStruct(estimatedGC_FN);
		trueGC.divideColBy(0,100);
		trueGC=averageMultipleGCs(trueGC);
		trueGC.write(trueGC_FN.replace(".txt", "_dupsAveraged.txt"));
		
		trueGC.keepRows(estimatedGC);
		trueGC = trueGC.mergeColumns(estimatedGC);
		trueGC.setColHeaders(new String[]{"true GC content","estimated GC content"});
		trueGC.write(writeFN);
	}

	private static MatrixStruct averageMultipleGCs(MatrixStruct trueGC) 
	{
		Hashtable<String, Double> gcHash = new Hashtable<String, Double>();
		for(int r = 0; r < trueGC.rows(); r++)
		{
			String rowName = trueGC.getRowHeaders()[r];
			if(!gcHash.containsKey(rowName))
			{
				gcHash.put(rowName, trueGC.matrix.get(r, 0));
			}
			else
			{
				double gcContent = (gcHash.get(rowName)+  trueGC.matrix.get(r, 0))/2;
				gcHash.put(rowName, gcContent);
			}	
		}
		
		MatrixStruct results = new MatrixStruct(gcHash.size(),1);
		Enumeration<String> sampleNames = gcHash.keys();
		int r = 0;
		while(sampleNames.hasMoreElements())
		{
			String sampleName = sampleNames.nextElement();
			results.setRow(r, sampleName, new double[]{gcHash.get(sampleName)});
			r++;
		}
		return results;
	}
}
