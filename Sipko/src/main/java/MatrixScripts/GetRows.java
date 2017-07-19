package MatrixScripts;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;

import Tools.FileUtils;
import Tools.Script;

public class GetRows extends Script<GetRows>
{
	String fileName = "E:/Groningen/Data/Juha/Genes31995/Healthy/PCA/31.07.pc1.illumina.genes.expressed_DownSamples/PC_1-0__zScores.txt";
	String fileName2 = "E:/Groningen/Data/Juha/Genes31995/Old/31.07.pc1.illumina.genes.expressed.DEseqnorm_notRounded/18DownSyndrome26Normal2Cancer_counts_transposed/PC_1-300_DevidedBySTdevsTop12000Expressed.txt";
	//String fileName2 = "ENSG00000268903,ENSG00000269981,ENSG00000225630";
	String writeName = "E:/Groningen/Data/Juha/Genes31995/Healthy/PCA/31.07.pc1.illumina.genes.expressed_DownSamples/PC_1-0_zScores_12000highest.txt";
	String remove = null;
	
	public void run()
	{
		try
		{
			MyMatrix file2 = null;
			
			if(!fileName2.contains(",") && fileName2.contains(".txt"))
				file2 = new MyMatrix(fileName2);
			else{
				String[] rowNames = fileName2.split(",");
				file2 = new MyMatrix(rowNames.length, 1);
				file2.rowNames = rowNames;
				file2.colNames[0] = "-";
			}
			
			BufferedReader reader = FileUtils.createReader(fileName);
			BufferedWriter writer = FileUtils.createWriter(writeName);
			String line = reader.readLine();
			writer.write(line+"\n");//write the first line by default (headers)
			Hashtable<String, Integer> toGet = file2.namesToHash(file2.rowNames);
			String[] results = new String[toGet.size()];
			while((line = reader.readLine()) != null)
			{	
				if(remove != null)
					line = line.replace(remove, "");
				String rowName = line.split("\t")[0];
				if(toGet.containsKey(rowName))
				{
					results[toGet.get(rowName)] = line;
					//writer.write(line+"\n");
				}
			}
			for(int r = 0; r < results.length; r++)
			{
				if(results[r] == null)
					continue;
				//System.out.println("r= "+(results[r]));
				writer.write(results[r]+"\n");
			}
			
			writer.close();
			reader.close();
			System.out.println("File written to:" + writeName);
		}catch(Exception e){e.printStackTrace();}
	}
}
