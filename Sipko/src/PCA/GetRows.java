package PCA;

public class GetRows 
{

	public static void main(String args[])
	{
		String fileName = "E:/Groningen/Data/PublicSamples/NoCancerSamples/est_counts_22214_samples.txt";
		String fileName2 = "E:/Groningen/Data/PublicSamples/DownSyndrome/18DownSyndrome26Normal274Cancer.txt";
		String writeName = "E:/Groningen/Data/PublicSamples/DownSyndrome/18DownSyndrome26Normal274Cancer_counts.txt";
		
		checkArgs(args);
		if(!System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
		{
			fileName = args[0];
			fileName2 = args[1];
			writeName = args[2];
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
		Matrix file1 = new Matrix(fileName);
		file1.keepRows(file2);
		file1.write(writeName);
		System.out.println("File written to:" + writeName);
	}
	public static void checkArgs(String[] args)
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
			return;
		if(args.length != 3)
		{
			System.out.println("Arguments supplied =" + args.length);
			System.out.println("This script retrieves rows for which rowNames are present in 1 file "
					+ "from another matrix file (automagically keeps header row).\n"
					+ "It uses the following 2 arguments:\n"
					+ "1. File to retriever rows from\n"
					+ "2.1 File with the rows to keep in the first column(header row is ignored)\n"
					+ "2.2 Alternatively you can use a comma separated list of rowNames you wish to retrive (SRR001,SR002,...)\n"
					+ "3. Name of the file to write to \n"
					+ "Make sure each file has at least 2 columns and rows (just at a bunch of 0's in the 2nd column if you must) \n");
			System.exit(1);
		}
	}
}
