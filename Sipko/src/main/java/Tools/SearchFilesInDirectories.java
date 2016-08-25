package Tools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import PCA.MatrixStruct;

public class SearchFilesInDirectories 
{
	//This script searches a folder (and all subdirectories) all for all files containing a certain "String" and
	//outputs all filenames (including pathname) that contain this string into a .txt file
	//same can be achieved with shell command : find . | grep "whateverYouAreSearchingFor" | grep -v "whateverYouWantToExclude"
	
	public static void main (String[] args) throws IOException
	{
		String folderName = "/groups/umcg-wijmenga/prm02/data_projects/Celiac_Disease/gene_expression/RNASEQ/RNAseq_NG_data_08032016/HiSeq/stimulated gluten specific Tcell clones";
		String writeName = "/local/groups/umcg-wijmenga/scr01/umcg-svandam/allSamples.txt";
		String[] searchString = new String[]{".fq"};
		String[] forbiddenString = new String[]{".md5"};
		
		//java -jar -Xmx1g SearchFilesInDirectories.jar folderName=/groups/umcg-pub/tmp04/public-rna-seq/ searchString=.fastq.gz writefn=/local/groups/umcg-bios/scr01/Sipko/Juha/samples.txt
		if(args.length == 0)
			checkArgs(args);
		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase()){
				case "foldername":
					folderName =value;
					break;
				case "searchstrings":
					searchString = value.split(",");
					break;
				case "writefn":
					writeName = value;
					break;
				case "forbiddenstrings":
					forbiddenString = value.split(",");
					break;	
				default:
					System.out.println("Incorrect argument supplied: "+ args[a] + "\n"
							+ "exiting");
					checkArgs(args);
					System.exit(1);
			}
		}
		if(args.length == 2)
			writeName = folderName+ searchString +".txt";
			
//		String folderName = "/groups/umcg-wijmenga/prm02/data_projects/Celiac_Disease/gene_expression/RNASEQ/RNAseq_NG_data_08032016/HiSeq/stimulated gluten specific Tcell clones";
//		String writeName = "/groups/umcg-wijmenga/prm02/data_projects/Celiac_Disease/gene_expression/RNASEQ/RNAseq_NG_data_08032016/HiSeq/stimulated gluten specific Tcell clones/allSamples.txt";
		
		
		File folder = new File(folderName);
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(writeName)));
		searchDirectory(folder, writer, searchString, forbiddenString);
		writer.close();
		System.out.println("done, file writen to:" + writeName);
	}
	public static void searchDirectory(File directory, String writeName, String[] searchString, String[] forbiddenString) throws IOException
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(writeName)));
		File[] files = directory.listFiles();
		
		for(File file : files)
		{
			if (file.isDirectory())
				searchDirectory(file, writer, searchString, forbiddenString);
			else 
			{
				boolean include = false;
				for(int s = 0; s < searchString.length; s++)
					if(file.getName().contains(searchString[s]))
						include = true;
				for(int s = 0; s < searchString.length; s++)
					if(file.getName().contains(forbiddenString[s]))
						include = false;
				if(include)
					writer.write(file.getAbsolutePath()+"\n");
			}
		}
		writer.close();
	}
	public static void searchDirectory(File directory, BufferedWriter writer, String[] searchString, String[] forbiddenString) throws IOException
	{
		File[] files = directory.listFiles();
		
		for(File file : files)
		{
			if (file.isDirectory())
				searchDirectory(file, writer, searchString, forbiddenString);
			else 
			{
				boolean include = false;
				for(int s = 0; s < searchString.length; s++)
					if(file.getName().contains(searchString[s]))
						include = true;
				for(int s = 0; s < forbiddenString.length; s++)
					if(file.getName().contains(forbiddenString[s]))
						include = false;
				
				if(include)
					writer.write(file.getAbsolutePath()+"\n");
			}
		}
	}
	public static ArrayList<String> searchDirectory(File directory, String searchString, String forbiddenString) throws IOException
	{
		File[] files = directory.listFiles();
		ArrayList<String> fastqFiles = new ArrayList<String>();
		for(File file : files)
		{
			if (file.isDirectory())
				searchDirectory(file, searchString, forbiddenString);
			else 
			{
				if(file.getName().contains(searchString) && !file.getName().contains(forbiddenString))
				{
					fastqFiles.add(file.getAbsolutePath());
				}
			}
		}
		return fastqFiles;
	}
	static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
			return;
		System.out.println("This script takes the following arguments:\n"
				+ "1. folderName=<foldername>, name of the folder to search through\n"
				+ "2. searchStrings=<searchString,searchString2>, a comma separated list of strings for which to include the files\n"
				+ "3. forbiddenStrings=<forbiddenString,forbiddenString2>, a comma separated list of strings that are not allowed to be in the filenames (default=.md5)\n"
				+ "4. writeFN=<writeFN>, name of the file to write (defaul=<searchString>.txt)\n");
		System.exit(1);
	}
}
