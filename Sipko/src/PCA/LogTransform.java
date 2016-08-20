package PCA;

import java.io.IOException;

public class LogTransform 
{
	//log 2 Transforms the data
	
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
		log2(writeFolder, expressionStruct, writeAll, 1);
	}

	public static void log2(String writeFolder, MatrixStruct expressionStruct, boolean writeAll, double add) throws IOException {
		JuhaPCA.PCA.log(" 9. Log2 transforming");
		expressionStruct.log2Transform(add);
		String quantFN = writeFolder+ "SAMPLE_NormalizedLog2.txt";
		
		JuhaPCA.PCA.log("10. Writing logged SAMPLE_Log2 normalized data in: " + quantFN);
		if(writeAll)expressionStruct.write(quantFN);
	}
	static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
			return;
		System.out.println("Wrong arguments");
		System.exit(1);
	}
		
}
