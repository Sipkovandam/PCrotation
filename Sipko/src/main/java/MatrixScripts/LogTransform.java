package MatrixScripts;

import java.io.File;
import java.io.IOException;

public class LogTransform 
{
	//log 2 Transforms the data
	static String filename = "";
	static String writeFolder = null;
	
	public static void main(String[] args) throws IOException 
	{
		boolean writeAll = true;
		checkArgs(args);
		if(writeFolder==null)
			writeFolder = new File(filename).getParent();
		MyMatrix expressionStruct = new MyMatrix(filename);
		log2(writeFolder, expressionStruct, writeAll, 1);
		String quantFN = writeFolder+ "/SAMPLE_NormalizedLog2.txt";
		System.out.println("File written to:" + quantFN);
	}

	public static void log2(String writeFolder, MyMatrix expressionStruct, boolean writeAll, double add) throws IOException {
		JuhaPCA.PCA.log(" 9. Log2 transforming");
		expressionStruct.log2Transform(add);
		String quantFN = writeFolder+ "SAMPLE_NormalizedLog2.txt";
		
		JuhaPCA.PCA.log("10. Writing logged SAMPLE_Log2 normalized data in: " + quantFN);
		if(writeAll)expressionStruct.write(quantFN);
	}
	static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko") && args.length < 1)
			return;
		if(args.length < 1)
		{
			System.out.println("Script requires the following argumetns:\n"
					+ "1. filename=<filename.txt> - File for which the values should be logged (log2)\n"
					+ "2. writefolder=<writefolder> - Name of path where output should be written (default=parentDirectory(<filename>)\n");
			System.exit(1);
		}
		
		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase())
			{
//				var = new JSONutil<Vars>().read(var.JSON_FN, var);
				case "filename":
					filename =value;
					break;
				case "writefolder":
					writeFolder =value;
					break;
				default:
					System.out.println("Incorrect argument supplied:\n"+ args[a] +"\nexiting");
					System.exit(1);
			}
		}
	}
		
}
