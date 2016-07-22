package PCA;

import java.io.IOException;

public class Spearman 
{
	//does not actually calculate the correlation, only ranks the sample
	public static void main(String[] args) throws IOException 
	{
		String expressionFN = "";
		String writeFolder = "";
		boolean writeAll = true;
		double spearman = 0;
		
		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase()){
				case "filename":
					expressionFN =value;
					break;
				case "writefolder":
					writeFolder = value;
					break;
				case "writeall":
					writeAll = Boolean.parseBoolean(value);
					break;
				default:
					checkArgs(args);
					System.out.println("Incorrect argument supplied; exiting");
					System.exit(1);
			}
		}
		
		MatrixStruct expressionStruct = new MatrixStruct(expressionFN);
		ranks(writeFolder, expressionStruct, spearman);
	}

	public static void ranks(String writeFolder, MatrixStruct expressionStruct, double spearman) throws IOException {
		String beforeRanks = writeFolder+"beforeRanks.txt";
		expressionStruct.write(beforeRanks);
		
		expressionStruct.expressionToRank(null,spearman);
		expressionStruct.write(writeFolder + "ranks.txt");
	}
	static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
			return;
		System.out.println("Wrong arguments");
		System.exit(1);
	}
}
