package PrepData;

import java.io.File;
import java.io.IOException;

public class KallistoAllSteps 
{
	public static String folder = null;
	public static String threshold = null;
	public static String transcriptsToGenesFile = null;
	public static String column = null;
	
	public static void main(String[] args) throws Exception
	{
		checkArgs(args);
		String fastqFilesFN = folder+"fastqs.txt";
		SearchFilesInDirectories.searchDirectory(new File(folder), fastqFilesFN, ".fastq", ".out");
		//RunKallistoOnFastqFiles();
		
	}
	static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
			return;
		if(args.length < 1)
		{
			System.out.println("Script requires the following argumetns:\n"
					+ "1. filename=<expressionFN.txt> - Expression file with genes on rows samples on columns\n"
					+ "2. writeFolder=<writeFolderFN.txt> - Folder where the files will be written (default=parentFolder(input.txt))\n"
					+ "3. geoFN=<geoFn.txt> - Optional file with geometric mean per gene to use (default=null)\n");
			System.exit(1);
		}
		
		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase())
			{
				case "folder":
					folder =value;
					break;
				case "threshold":
					threshold = value;
					break;
				case "tstogenefn":
					transcriptsToGenesFile = value;
					break;
				case "column":
					column = value;
					break;
				default:
					System.out.println("Incorrect argument supplied:\n"+ args[a] +"\nexiting");
					System.exit(1);
			}
		}
	}
}
