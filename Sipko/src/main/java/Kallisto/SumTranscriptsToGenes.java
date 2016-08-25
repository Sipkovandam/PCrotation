package Kallisto;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.Hashtable;

import JuhaPCA.PCA;
import PCA.Matrix;
import Tools.FileUtils;

public class SumTranscriptsToGenes 
{
	//Sums up all read counts originating from the same gene
	//transcripts should be on rows, but script will transpose if needed
	static String countFileName = "E:/Groningen/Data/LifeLines/Counts.txt";
	static String transcriptToGeneFile = "E:/Groningen/Data/Annotation/hg19.v75.cdna.all.enst2ensg.txt";
	static String writeName = null;
	
	public static void main(String[] args)
	{
		checkArgs(args);
		
		run(countFileName, transcriptToGeneFile);
		
		System.out.println("Done!");
	}

	static Matrix run(String countFileName, String transcriptToGeneFile) 
	{
		System.out.println(countFileName);
		if(writeName == null)
			writeName = FileUtils.replaceEnd(countFileName, "_GENES.txt.gz");
		
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

	static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko") && args.length < 1)
			return;
		if(args.length < 1)
		{
			System.out.println("Script requires the following argumetns:\n"
					+ "1. countFN=<countfilename.txt> - File with the counts for every sample\n"
					+ "2. transcriptToGeneFN=<transcripttogenefile.txt> - File with transcript names in first column and gene names in second\n"
					+ "3. writeFN=<writename.txt> - Name of the file to be written (default=<countFN>+_GENES.txt.gz\n");
			System.exit(1);
		}
		
		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase())
			{
//				var = new JSONutil<Vars>().read(var.JSON_FN, var);
				case "countfn":
					countFileName =value;
					break;
				case "transcripttogenefn":
					transcriptToGeneFile =value;
					break;
				case "writefn":
					writeName =value;
					break;
				default:
					System.out.println("Incorrect argument supplied:\n"+ args[a] +"\nexiting");
					System.exit(1);
			}
		}
	}
}
