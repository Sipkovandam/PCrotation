package PCA;

import java.io.IOException;

public class QuantileNormalize 
{
	public static void main(String[] args) throws IOException 
	{
		String expressionFN = "";
		String writeFolder = "";
		boolean writeAll = true;

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
		quantileNormalize(expressionStruct, writeFolder, writeAll);
	}

	public static void quantileNormalize(MatrixStruct expressionStruct, String writeFolder, boolean writeAll) throws IOException {
		pca.PCA.log(" 6. Calculating quantile normalization vector");
		MatrixStruct qNormVector = expressionStruct.quantileNormVector();
		qNormVector.write(writeFolder+ "SAMPLE_QuantileVector.txt");
	
		pca.PCA.log(" 7. Quantile normalization");
		expressionStruct.expressionToRank(qNormVector,0);
		String quantFNnotLogged = writeFolder+ "SAMPLE_QuantileNormalized.txt";
		pca.PCA.log(" 8. Writing quantile normalized data in: " + quantFNnotLogged);
		if(writeAll)expressionStruct.write(quantFNnotLogged);
	}
	static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
			return;
		System.out.println("Wrong arguments");
		System.exit(1);
	}
}
