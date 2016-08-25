package PrepData;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import PCA.MatrixStruct;

public class GeneToTranscript 
{

	public static void main(String[] args) throws IOException
	{
		String genesFN = "E:/Groningen/Data/GenePositionInfo_23X_24Y_25MT_26rest.txt";
		String transcriptsToGenesFN = "E:/Groningen/Data/Annotation/hg19.v75.cdna.all.enst2ensg.txt";
		String writeFN = genesFN.replace("Gene", "Transcript");
		
		BufferedReader reader = new BufferedReader(new FileReader(new File(transcriptsToGenesFN)));
		String line = null;
		Hashtable<String,String> transcriptsToGenes = new Hashtable<String,String>();
		ArrayList<String> transcripts = new ArrayList<String>();
		while((line = reader.readLine()) != null)
		{
			String[] eles = line.split("\t");
			transcriptsToGenes.put(eles[0], eles[1]);
			transcripts.add(eles[0]);
		}
		reader.close();
		
		MatrixStruct genes = new MatrixStruct(genesFN);
		
		MatrixStruct output = new MatrixStruct(transcriptsToGenes.size(),3);
		output.setColHeaders(new String[]{"Chromosome","Start","End"});
		output.setRowHeaders(transcripts.toArray(new String[transcripts.size()]));
		transcriptsToGenes.keys();
		
		int adj = 0;
		for(int r = 0; r < output.rows(); r++)
		{
			String transcript = output.getRowHeaders()[r];
			String geneName = transcriptsToGenes.get(transcript);
			System.out.println("geneName = " + geneName + " " + transcriptsToGenes.size() + " " + transcript);
			if(!genes.rowHash.containsKey(geneName))
			{
				adj++;
				continue;
			}
			int row = genes.rowHash.get(geneName);
			double chr = genes.matrix.get(row, 0);
			double start = genes.matrix.get(row, 1);
			double end = genes.matrix.get(row, 2);
			output.setRow(r-adj, transcript, new double[]{chr,start,end});
		}
		output.write(writeFN);
		System.out.println("File written to:" + writeFN);
	}
}

