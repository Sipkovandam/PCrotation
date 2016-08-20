package PCA;

public class Transpose 
{
	//Transposes the input file.
	
	public static void main (String[] args)
	{
		//String fileName = "E:/Groningen/Data/PublicSamples/05-2016/22214Samples_Voom_Correl/SAMPLE_NormalizedLog2_10Rows.txt";//
		String fileName = "E:/Groningen/Data/Juha/Rlog67000Samples/31.07.expressed100Lines.txt";
		if(args.length < 1)
			checkArgs(args);
		if(!System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
		{
			for(int a = 0; a < args.length; a++)
			{
				String arg = args[a].split("=")[0];
				String value = args[a].split("=")[1];
				switch (arg.toLowerCase()){
					case "filename":
						fileName = value;
						break;
					default:
						checkArgs(args);
						System.out.println("Incorrect argument supplied; exiting");
						System.exit(1);
				}
			}
		}
		
		Matrix matrix = new Matrix(fileName);
		//matrix.print(5,5);
		System.out.println("Rows = " + matrix.rowNames.length + " columns = " + matrix.colNames.length);
		matrix.writeTransposed(fileName.replace(".txt", "_transposed.txt"));
		
//		matrix.write(fileName.replace(".txt", "_first10.txt"), -1,10);
//		matrix.write(fileName.replace(".txt", "_first100geneEvs.txt"), -1,100);
//		matrix.write(fileName.replace(".txt", "_transposed_100.txt"), 100, -1);
//		matrix.write(fileName.replace(".txt", "_transposed_300.txt"), 300, -1);
//		matrix.write(fileName.replace(".txt", "_transposed_500.txt"), 500, -1);
//		matrix.write(fileName.replace(".txt", "_transposed_1000.txt"), 1000, -1);
		System.out.println("Done");
	}
	
	public static void checkArgs(String[] args)
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
			return;
		if(args.length < 1)
		{
			System.out.println("Arguments supplied =" + args.length);
			System.out.println("This script transposes a file and needs the following argument:\n"
					+ "1. fileName=<fileName> -  File to retrieve rows from\n"
					+ "Example: java -jar -Xmx32G filename=filename.txt\n");
			System.exit(1);
		}
	}
}
