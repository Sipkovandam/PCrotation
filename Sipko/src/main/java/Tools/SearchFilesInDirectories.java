package Tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import PCA.MatrixStruct;

public class SearchFilesInDirectories 
{
	//This script searches a folder (and all subdirectories) all for all files containing a certain "String" and
	//outputs all filenames (including pathname) that contain this string into a .txt file
	//same can be achieved with shell command : find . | grep "whateverYouAreSearchingFor" | grep -v "whateverYouWantToExclude"
	
	public static void main (String[] args) throws IOException
	{	
		String folders = null;//"E:/Groningen/Test/STAR/STAR/,E:/Groningen/Test/PCA/Test2/";
		String writeName = null;//"E:/Groningen/Test/STAR/STAR/textsearch.txt";
		String[] searchString = new String[]{".txt"};
		String requiredStringFN = null;//"E:/Groningen/Test/STAR/STAR/RequiredStrings.txt";
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
					folders =parseString(value);
					break;
				case "folders":
					folders =parseString(value);
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
				case "requiredstringfn":
					requiredStringFN = parseString(value);
					break;	
				default:
					System.out.println("Incorrect argument supplied: "+ args[a] + "\n"
							+ "exiting");
					checkArgs(args);
					System.exit(1);
			}
		}
		if(args.length == 2 || writeName==null || writeName.length()<2)
			writeName = folders.split(",")[0]+ searchString +".txt";
			
//		String folderName = "/groups/umcg-wijmenga/prm02/data_projects/Celiac_Disease/gene_expression/RNASEQ/RNAseq_NG_data_08032016/HiSeq/stimulated gluten specific Tcell clones";
//		String writeName = "/groups/umcg-wijmenga/prm02/data_projects/Celiac_Disease/gene_expression/RNASEQ/RNAseq_NG_data_08032016/HiSeq/stimulated gluten specific Tcell clones/allSamples.txt";
		
		String[] folderNames = folders.split(",");
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(writeName)));
		final String[] searchStringFinal = searchString;
		final String[] forbiddenStringFinal = forbiddenString;
		final String[] requiredStringFinal = readRequiredStrings(requiredStringFN);
		Stream.of(folderNames).forEach(folderName -> searchDirectory(new File(folderName), writer, searchStringFinal,requiredStringFinal, forbiddenStringFinal));
		writer.close();
		System.out.println("done, file writen to:" + writeName);
	}
	private static String[] readRequiredStrings(String requiredStringFN) throws FileNotFoundException, IOException 
	{
		if(requiredStringFN == null)
			return null;
		BufferedReader reader = FileUtils.createReader(requiredStringFN);
		String[] requiredStrings = reader.lines().map(line -> 
		{
			return line;
		}).collect(Collectors.toList()).toArray(new String[]{});

		return requiredStrings;
	}
	public static void searchDirectory(File directory, BufferedWriter writer, String[] searchString, String[] forbiddenString)
	{
		searchDirectory(directory, writer, searchString, null, forbiddenString);
	}
	public static void searchDirectory(File directory, BufferedWriter writer, String[] searchString, String[] requiredString, String[] forbiddenString)
	{
		
		try 
		{
			File[] files = directory.listFiles();
			for(File file : files)
			{
				if (file.isDirectory())
					searchDirectory(file, writer, searchString,requiredString, forbiddenString);
				else 
				{
					String fileName = file.getName();
					boolean include = checkInclude(searchString, fileName, requiredString, forbiddenString);
					
					if(include)
							writer.write(file.getAbsolutePath()+"\n");	
				}
			}
		} catch (IOException e) {e.printStackTrace();}
	}
	private static boolean checkInclude(String[] searchString, String fileName, String[] requiredString, String[] forbiddenString) {
		for(int s = 0; s < forbiddenString.length; s++)
			if(fileName.contains(forbiddenString[s]))
				return false;
		
		for(int s = 0; s < searchString.length; s++)
			if(fileName.contains(searchString[s]))
				if(requiredString == null)
					return true;
				else//check if the other required string (any of them) is also present in the filename
					for(int r = 0; r < requiredString.length; r++)
						if(fileName.contains(requiredString[r]))
							return true;
		return false;
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
	private static String parseString(String value) {
		if(value.equals("null"))
			return null;
		return value;
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
