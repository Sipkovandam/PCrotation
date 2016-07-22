package PCA;

import java.io.IOException;

public class CorrectReadcounts {

	public static void main(String[] args) throws IOException 
	{
		String expressionFN = "";
		String writeFolder = "";
		boolean writeAll = true;	
		double correctTotalReadCount = 1000000;
		double add = 0.5;

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
				case "correcttotalreadcount":
					correctTotalReadCount = Double.parseDouble(value);
					break;	
				case "add":
					add = Double.parseDouble(value);
					break;	
				default:
					checkArgs(args);
					System.out.println("Incorrect argument supplied; exiting");
					System.exit(1);
			}
		}
		
		MatrixStruct expressionStruct = new MatrixStruct(expressionFN);
		correct(writeFolder, correctTotalReadCount, expressionStruct, writeAll, add);
	}

	public static void correct(String writeFolder, double correctTotalReadCount, MatrixStruct expressionStruct, boolean writeAll, double add) throws IOException {
		JuhaPCA.PCA.log(" 6. Correcting for total read count");
		String correctedNotLogged =  writeFolder+ "SAMPLE_TotalReadCountNormalized.txt";
		expressionStruct.correctForTotalReadCount(correctTotalReadCount, add);
		if(writeAll)expressionStruct.write(correctedNotLogged);
	}
	static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
			return;
		System.out.println("Wrong arguments");
		System.exit(1);
	}
}
