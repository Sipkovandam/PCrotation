package PCA;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Hashtable;

import MatrixScripts.MyMatrix;

public class SortChromosome 
{
	//Sorts genes in a file by chromosome location.
	
	public static void main(String[] args) throws IOException
	{
		String sampleFile = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/RandomSamples_Adj/" + "Adjusted_PC0_transposed.txt";
		String chromLocationsFile = "E:/Groningen/Data/GenePositionInfo_23X_24Y_25MT_26rest.txt";
		boolean writeAll = true;
		String writeFolder = null;
		
		checkArgs(args);
		if(!System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko") && args.length !=0)
		{
			sampleFile = args[0];
			chromLocationsFile = args[1];
		}
		sortChromosome(sampleFile, chromLocationsFile, writeAll, writeFolder);
		
	}
	public static void checkArgs(String[] args)
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
			return;
		if(args.length!= 2)
		{
			System.out.println("This scripts sorts gene IDs in the rowNames based on genomic location:\n"
					+ "It uses the following 2 arguments:\n"
					+ "1. File with gene names on rows \n"
					+ "2. File with:"
					+ "		- 1st column: gene names\n"
					+ "   	- 2nd column: Chromosome (numbers only, replace X,Y etc with numbers first)\n"
					+ "		- 3rd column: Genomic location\n");
			System.exit(1);
		}
	}
	public static MyMatrix sortChromosome(String sampleFile, String chromLocationsFile, boolean writeAll, String writeFolder) throws IOException
	{
		MyMatrix sampleStruct = new MyMatrix(sampleFile);
		sampleStruct.putGenesOnRows();
		sampleStruct = SortChromosome.sort(sampleStruct , chromLocationsFile);
		sampleFile = writeFolder+"chromSorted.txt";
		if(writeAll)
		{
			sampleStruct.write(sampleFile);
			System.out.println("File written to: " + writeFolder+"chromSorted.txt");
		}
		return sampleStruct;
	}
	public static MyMatrix sort(MyMatrix samples, String chromLocationsFile)
	{
		if(chromLocationsFile==null)
			return samples;
		
		File chromfile = new File(chromLocationsFile);
		if(!chromfile.exists())
		{
			System.out.println("Genes are not sorted based on their chromosome locations because this file does not exist: \n" + chromfile.getAbsolutePath());
			return samples;
		}
	
		JuhaPCA.PCA.log("Sorting IDs based on chromosome locations");
		MyMatrix chromLocs = new MyMatrix(chromLocationsFile);
		//only sorts those that are in the sample
		samples.keepRows(chromLocs);
		samples.rows();
		System.out.println("IDs present in both files = " + chromLocs.rows());
		//first sort the chromosome location file
		String[] sortedNames = sortChromMatrix(chromLocs);
		Hashtable<String, Integer> geneIndex = MyMatrix.makeHash(sortedNames);
		MyMatrix sortedSample = new MyMatrix(samples.rows(),samples.cols());
		sortedSample.setColHeaders(samples.getColHeaders());
		int unsorted = 0;
		for(int r = 0; r < samples.rows(); r++)
		{ 
			int newRowNumber = 0;
			
			if(geneIndex.get(samples.getRowHeaders()[r]) != null)
			{
				newRowNumber = geneIndex.get(samples.getRowHeaders()[r]);
			}
			else//if it is not in the file that determines the order just put it at the end;
			{
				System.out.println("This ID is not present in the sort file" + samples.getRowHeaders()[r]);
				newRowNumber = samples.getRowHeaders().length-unsorted-1;
				unsorted++;
			}
			//System.out.println(sample.getRowHeaders()[r].length());
			sortedSample.setRow(newRowNumber,samples.getRowHeaders()[r], samples.getRowValues(r));
		}
		System.out.println("Missing " + unsorted + " IDs in the file used for the sort order."
				+ "The rows for which IDs were missing have been added at the end of the sorted list");
		
		return sortedSample;
	}
	private static String[] sortChromMatrix(MyMatrix chromLocs) 
	{
		class Row
		{
			public Row(String geneName1, double location1) {
				geneName = geneName1;
				location = location1;
			}
			String geneName;
			double location;
		};
		Hashtable<Integer, ArrayList<Row>> genome = new Hashtable<Integer, ArrayList<Row>>();
		
		int lastChrom = 0;
		//put all the genes in arraylists according to their chromosome number
		for(int r = 0; r < chromLocs.rows(); r++)
		{
			int chromNumber = (int) chromLocs.matrix.get(r, 0);
			ArrayList<Row> chrom = genome.get(chromNumber);
			if(chrom == null)
				chrom = new ArrayList<Row>();
			chrom.add(new Row(chromLocs.getRowHeaders()[r], chromLocs.matrix.get(r, 1)));
			genome.put(chromNumber, chrom);
				
			if(lastChrom < chromNumber)
				lastChrom = chromNumber;
		}
		int r = 0;
		String[] orderedNames = new String[chromLocs.rows()];
		for(int c = 1; c < lastChrom+1; c++)
		{
			ArrayList<Row> chrom = genome.get(c);
			if(chrom == null)//can happen when all genes on a chromosome are not present in the target file
				continue;
			//sort on chromosome location (start pos)
			Collections.sort(chrom, new Comparator<Row>()
			{
				@Override
				public int compare(Row s1, Row s2)
				{
					if (s1.location > s2.location)
						return 1;	// tells Arrays.sort() that s1 comes after s2
					else if (s1.location < s2.location)
						return -1;   // tells Arrays.sort() that s1 comes before s2
					else 
						return 0;
				}
			});
			for(int gene = 0; gene < chrom.size(); gene++)
			{
				orderedNames[r] = chrom.get(gene).geneName;
				r++;
			}
		
		}
		return orderedNames;
	}
}
