package PrepData;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Hashtable;
import java.io.FileWriter;

import pca.MatrixStruct;

public class RemoveBadSamples {
	//removes samples with expression for less then 10% of the genes
	
	public static void main (String[] args) throws IOException
	{
		String fn = args[0];
		String newName = "_removedBadSamples.txt.gz";
		String writeFN = fn.replace(".txt", "").replace(".gz", "") + newName;
		MatrixStruct expression = new MatrixStruct(fn);
		int cutOff = expression.rows()/10;
		
		Hashtable<String,Integer> genesExpressed = new Hashtable<String,Integer>();
		expression.putGenesOnRows();
		int nGood = 0;
		
		//determine which samples are good
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(fn.replace(".txt", "").replace(".gz", "") + "genesExpressed.txt")));
		writer.write("Sample\tExpressed\tTotal");
		for(int c = 0; c < expression.cols(); c++)
		{
			int expressed = 0;
			for(int r = 0; r < expression.rows(); r++)
			{
				if(expression.matrix.get(r,c) > 0)
					expressed++;
			}
			writer.write(expression.getColHeaders()[c]+"\t"+expressed+"\t"+expression.rows()+"\n");
			genesExpressed.put(expression.getColHeaders()[c], expressed);
			if(expressed > cutOff)
				nGood++;
			else
				System.out.println(expression.getColHeaders()[c]);//samples that do not make the cutoff (e.g. the bad samples)
		}
		writer.close();
		MatrixStruct result = new MatrixStruct(expression.rows(),nGood);
		result.setRowHeaders(expression.getRowHeaders());
		
		int outCol = 0;
		for(int c = 0; c < expression.cols(); c++)
		{
			if(genesExpressed.get(expression.getColHeaders()[c]) <= cutOff)
				continue;
			result.setColHeader(outCol,expression.getColHeaders()[c]);
			for(int r = 0; r < result.rows(); r++)
			{
				result.matrix.set(r, outCol, expression.matrix.get(r, c));
			}
			outCol++;
		}
		result.write(writeFN);
		System.out.println("Done, file written to:" + writeFN);
	}
}
