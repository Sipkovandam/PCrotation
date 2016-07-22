package PrepData;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.Hashtable;

import PCA.Matrix;
import pca.PCA;

public class SumTranscriptsToGenes 
{
	//Sums up all read counts originating from the same gene
	//transcripts should be on rows, but script will transpose if needed
	
	public static void main(String[] args)
	{
		String countFileName = "E:/Groningen/Data/LifeLines/Counts.txt";
		String transcriptToGeneFile = "E:/Groningen/Data/Annotation/hg19.v75.cdna.all.enst2ensg.txt";
		
		if(args.length == 2)
		{
			countFileName = args[0];
			transcriptToGeneFile = args[1];
			//System.out.println("2 arguments passed to function");			
		}
		
		run(countFileName, transcriptToGeneFile);
		
		System.out.println("Done!");
	}

	static Matrix run(String countFileName, String transcriptToGeneFile) 
	{
		System.out.println(countFileName);
		String writeName = countFileName.replace(".txt", "GENES.txt");
		Hashtable<String,String> transcriptToGeneIDs = null;
		Hashtable<String,Integer> transcriptsPerGeneCounts = null;
		
		// create the conversion hashtable
		try 
		{
			BufferedReader reader = new BufferedReader(new FileReader(new File(transcriptToGeneFile)));
			String line = null;
			transcriptToGeneIDs = new Hashtable<String,String>();
			transcriptsPerGeneCounts = new Hashtable<String,Integer>(); //counts how many different transcript IDs there for each geneID
			
			while((line = reader.readLine())!=null)
			{
				String[] eles = line.split("\t");
				transcriptToGeneIDs.put(eles[0], eles[1]);
				
				//Count how many transcripts map to each gene
				if(transcriptsPerGeneCounts.containsKey(eles[1]))
					transcriptsPerGeneCounts.put(eles[1], transcriptsPerGeneCounts.get(eles[1])+1);
				else
				{
					transcriptsPerGeneCounts.put(eles[1], 1);
				}
			}
			

			reader.close();
		} catch (FileNotFoundException e) {e.printStackTrace();} catch (IOException e) {e.printStackTrace();}
		
		
		Matrix counts = new Matrix(countFileName);//matrix with the counts of the transcripts		Matrix counts = new Matrix(countFileName);//matrix with the counts of the transcripts
		System.out.println("First row name:" + counts.rowNames[0] + "First column name:" + counts.colNames[0]);
		if(counts.colNames[0].contains("ENSG") || counts.colNames[0].contains("ENST"))
		{
			System.out.println("RowNames do not contain ensembl IDs in original matrix; Transposing");
			counts.transpose();
		}
		//create a hashtable that has the geneNames and the row of these names in the results file
		Hashtable<String,Integer> geneIDtoRow = new Hashtable<String,Integer>();
		
		String[] rowNames = new String[transcriptsPerGeneCounts.size()];
		Enumeration<String> keys = transcriptsPerGeneCounts.keys();
		System.out.println("There are "+ transcriptsPerGeneCounts.size() + " different gene IDs in this file");
		String key = null;
		int r = 0;
		while(keys.hasMoreElements())
		{
			key = keys.nextElement();
			rowNames[r] = key;
			r++;
		}
		
		Arrays.sort(rowNames);
		System.out.println(rowNames.length+","+counts.colNames.length);
		Matrix results = new Matrix(rowNames.length, counts.colNames.length);
		System.out.println("rows =" + rowNames.length + " cols =" + counts.colNames.length);
		results.rowNames = rowNames;
		results.colNames = counts.colNames;
		
		for(int x = 0; x < results.rowNames.length; x++)
		{
			geneIDtoRow.put(results.rowNames[x], x);
		}
		
		//Count reads per gene for the genes you want to keep
		System.out.println("Time halfway point:");
		
		Matrix averageExpression = new Matrix(rowNames.length, 2);
		averageExpression.colNames[0] = "Average Expression";
		averageExpression.colNames[1] = "Number of different transcripts";
		averageExpression.rowNames = results.rowNames;
		
		//for each column do the summing
		for(int y = 0; y < counts.colNames.length;y++)
		{
			if(y%100==0)
			{
				PCA.log("Progress:" + y +"/"+ counts.colNames.length);
			}
			for(int x = 0; x < counts.rowNames.length;x++)//add the counts of the transcripts to each of the corresponding genes for this column
			{
				String geneID = transcriptToGeneIDs.get(counts.rowNames[x]);
				if(geneID== null)
					continue;
				double n = counts.values[x][y];
				
				int row = 0;
				if(geneIDtoRow.containsKey(geneID))
					row = geneIDtoRow.get(geneID);
				else
					continue;
				results.values[row][y]+=n;
				averageExpression.values[row][0]+=n;
				averageExpression.values[row][1]++;
			}
		}
		for(int row = 0; row < averageExpression.rowNames.length; row++)
		{
			averageExpression.values[row][0]/=counts.colNames.length;
			averageExpression.values[row][1]/=counts.colNames.length;
		}
		averageExpression.write(writeName.replace(".txt", "_Averages.txt"));
		results.print(9,5);
		
		results.write(writeName, -1);//-1 indicates all decimals should be written
		return results;
	}

	private static void printMsg() {
		System.out.println("This function takes 2 arguments");
		
	}

}
