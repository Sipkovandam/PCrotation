package Analyses;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;

import MatrixScripts.MatrixStruct;

public class IrisSampleNamesAdjust 
{
	//changes colnames to the folder name the files were contained in, as this tells what kind of sample it is (treated/untreated,biopsy etc...)
	public static void main(String[] args) throws IOException
	{
		String sampleNamesFN = "E:/Groningen/Data/Iris/allSamples_edited.txt";
//		String samplesFN = "E:/Groningen/Data/Iris/CountsSD200/CountsGENESIrisColsEdited.txt";
//		String writeFN = "E:/Groningen/Data/Iris/CountsSD200/CountsGENESIrisCorrectNames.txt";
		
		String samplesFN = "E:/Groningen/Data/Iris/Merged/CountsGENES_Merged.txt";
		String writeFN = "E:/Groningen/Data/Iris/Merged/CountsGENES_Merged_CorrectNames.txt";
		MatrixStruct samples = new MatrixStruct(samplesFN);
		
		BufferedReader reader = new BufferedReader(new FileReader(new File(sampleNamesFN)));
		Hashtable<String,String> oldToNew = new Hashtable<String,String>();
		String line = null;
		
		int x = 0;
		while((line= reader.readLine()) != null)
		{
			if(!line.contains("/"))
				continue;
			File file = new File(line);
			String currentName = file.getName().replace("_1.fq.gz", "").replace("_2.fq.gz", "").replace("_R1_", "").replace("_R2_", "").replace(".fq.gz", "");
			String parentName = file.getParent();
			String name = parentName.substring(parentName.lastIndexOf("\\")+1, parentName.length());
			if(line.toLowerCase().contains("discarded"))
				name = "DISCARDED_"+ name;
			System.out.println("name = " + currentName + " output name " + name);
			oldToNew.put(currentName,name);
		}
		
		for(int c = 0; c < samples.cols(); c++)
		{
			String colName = samples.getColHeaders()[c];
			System.out.println(colName);
			File file = new File(colName);
			colName = file.getName();
			System.out.println(colName);
			samples.setColHeader(c,oldToNew.get(colName));
		}
		reader.close();
		samples.write(writeFN);
		System.out.println("FileName written in " + writeFN);
	}
}
