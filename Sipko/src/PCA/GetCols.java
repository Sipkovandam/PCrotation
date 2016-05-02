package PCA;

import java.util.Hashtable;

public class GetCols {

	public static void main(String[] args)
	{
		String fileName = "E:/Groningen/Data/PublicSamples/Test13/directPCA_Voom_0.2/est_counts_nocancernocellline/PC_1-300_DevidedBySTdevs.txt";
		//String fileName2 = "E:/Groningen/Data/PublicSamples/Test9/PublicSamplesWithoutDownSyndrome.txt";
		String fileName2 = "SRR514129";
		String writeName = "E:/Groningen/Data/PublicSamples/Test13/directPCA_Voom_0.2/est_counts_nocancernocellline/PC_1-300_DevidedBySTdevs_SRR514129.txt";
		String remove = "_200";
		
		checkArgs(args);
		if(args.length!=0)
		{
			for(int a = 0; a < args.length; a++)
			{
				String arg = args[a].split("=")[0];
				String value = args[a].split("=")[1];
				switch (arg.toLowerCase()){
					case "filename":
						fileName = value;
						break;
					case "columnstoget":
						fileName2 = value;
						break;
					case "getgenes":
						fileName2 = value;
						break;
					case "writename":
						writeName = value;
						break;
					case "remove":
						remove = value;
						break;
					default:
						checkArgs(args);
						System.out.println("Incorrect argument supplied; exiting");
						System.exit(1);
				}
			}
		}
	
		Matrix file2 = null;

		if(!fileName2.contains(",") && fileName2.contains(".txt"))
			file2 = new Matrix(fileName2);
		else{
			String[] rowNames = fileName2.split(",");
			file2 = new Matrix(rowNames.length, 1);
			file2.rowNames = rowNames;
			file2.colNames[0] = "-";
		}
	
		Matrix samples = new Matrix(fileName);
		for(int c = 0; c < samples.colNames.length; c++)
			samples.colNames[c] = samples.colNames[c].replace(remove, "");
		
		Hashtable<String, Integer> sampleColHash = samples.colNamesToHash();
		
		int n = 0;
		for(int s = 0; s < file2.rowNames.length; s++)
		{
			if(sampleColHash.containsKey(file2.rowNames[s]))
			{
				n++;
			}
		}
		
		Matrix output = new Matrix(samples.rowNames.length, n);
		output.rowNames = samples.rowNames;
	
		int o = 0;
		for(int s = 0; s < file2.rowNames.length; s++)
		{
			if(sampleColHash.containsKey(file2.rowNames[s]))
			{
				int col = sampleColHash.get(file2.rowNames[s]);
				output.colNames[o] = file2.rowNames[s];
				
				for(int r = 0; r< samples.rowNames.length; r++)
					output.values[r][o] = samples.values[r][col];
				o++;
			}
		}
		output.write(writeName);
		System.out.println("Done! File written to:" + writeName);
	}
	public static void checkArgs(String[] args)
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
			return;
		if(args.length < 3)
		{
			System.out.println("Arguments supplied =" + args.length);
			System.out.println("This script retrieves columns for which rowNames are present in the 'getGenes file' "
					+ "from another matrix file (automagically keeps header row).\n"
					+ "It uses the following 2 arguments:\n"
					+ "1. fileName=<fileName> -  File to retrieve rows from\n"
					+ "2.1 columnsToGet=<fileName> -  File with the cols to keep in the first column(header row is ignored)\n"
					+ "2.2 columnsToGet=<gene1,gene2,gene3> - Alternatively you can use a comma separated list of rowNames you wish to retrive (SRR001,SR002,...)\n"
					+ "3. writeName=<fileName> - Name of the file to write to \n"
					+ "Make sure each file has at least 2 columns and rows (just at a bunch of 0's in the 2nd column if you must) \n");
			System.exit(1);
		}
	}
}
