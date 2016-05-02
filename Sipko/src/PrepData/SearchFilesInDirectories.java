package PrepData;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import pca.MatrixStruct;

public class SearchFilesInDirectories 
{

	public static void main (String[] args) throws IOException
	{
		String folderName = "/groups/umcg-wijmenga/prm02/data_projects/Celiac_Disease/gene_expression/RNASEQ/RNAseq_NG_data_08032016/HiSeq/stimulated gluten specific Tcell clones";
		String writeName = "/local/groups/umcg-wijmenga/scr01/umcg-svandam/allSamples.txt";
		String searchString = ".fq";
		String forbiddenString = ".md5";
		
		if(args.length >0)
		{
			folderName = args[0];
			searchString = args[1];
			writeName = folderName+ searchString +".txt";
		}
//		String folderName = "/groups/umcg-wijmenga/prm02/data_projects/Celiac_Disease/gene_expression/RNASEQ/RNAseq_NG_data_08032016/HiSeq/stimulated gluten specific Tcell clones";
//		String writeName = "/groups/umcg-wijmenga/prm02/data_projects/Celiac_Disease/gene_expression/RNASEQ/RNAseq_NG_data_08032016/HiSeq/stimulated gluten specific Tcell clones/allSamples.txt";
		
		
		File folder = new File(folderName);
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(writeName)));
		searchDirectory(folder, writer, searchString, forbiddenString);
		writer.close();
		System.out.println("done, file writen to:" + writeName);
	}
	public static void searchDirectory(File directory, BufferedWriter writer, String searchString, String forbiddenString) throws IOException
	{
		File[] files = directory.listFiles();
		
		for(File file : files)
		{
			if (file.isDirectory())
				searchDirectory(file, writer, searchString, forbiddenString);
			else 
			{
				if(file.getName().contains(searchString) && !file.getName().contains(forbiddenString))
				{
					writer.write(file.getAbsolutePath()+"\n");
				}
			}
		}

	}
}
