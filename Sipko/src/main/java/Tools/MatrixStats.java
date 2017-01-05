package Tools;

import java.io.IOException;

import PCA.Matrix;
import PCA.MatrixStruct;

public class MatrixStats {

	public static String matrixFN = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/TESTexpression.txt";
	public static String writeFN_Rows = null;
	public static String writeFN_Cols = null;
	
	public static void main(String[] args) throws IOException 
	{
		checkArgs(args);
		//rows
		MatrixStruct matrix = new MatrixStruct(matrixFN);
		MatrixStruct averagesRows = matrix.getAveragesPerRow();
		MatrixStruct stdevsRows = matrix.stDevRows();
		MatrixStruct mergeRows = averagesRows.mergeColumns(stdevsRows);
				
		if(writeFN_Rows == null)
			if(mergeRows.getRowHeaders()[0].contains("ENSG0") || mergeRows.getRowHeaders()[0].contains("ENST0"))
			{
				writeFN_Rows=FileUtils.replaceEnd(matrixFN, "_Gene_stats.txt");
				writeFN_Cols=FileUtils.replaceEnd(matrixFN, "_Sample_stats.txt");
			}
			else
			{
				writeFN_Rows=FileUtils.replaceEnd(matrixFN, "_Row_stats.txt");
				writeFN_Cols=FileUtils.replaceEnd(matrixFN, "_Col_stats.txt");
			}

		mergeRows.write(writeFN_Rows);
		
		//cols
		MatrixStruct colAverages = matrix.getAveragesPerCol();
		MatrixStruct colSstdevs = matrix.stDevCols();
		MatrixStruct colMerge = colAverages.mergeColumns(colSstdevs);
		colMerge.write(writeFN_Cols);
	}
	static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko") && args.length < 1)
			return;
		if(args.length < 1)
		{
			System.out.println("This scritps calculates the averages and stdevs per row and column\n"
					+ "Script requires the following arguments:\n"
					+ "1. matrix=<matrix.txt> - input matrix to get stats for\n"
					+ "2. writeFNRows=<writeFileName.txt> - optional\n"
					+ "3. writeFNCols=<writeFileName.txt> - optional\n");
			System.exit(1);
		}
		
		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase())
			{
				case "matrix":
					matrixFN =value;
					break;
				case "writefnrows":
					writeFN_Rows =value;
					break;
				case "writefncols":
					writeFN_Rows =value;
					break;
				default:
					System.out.println("Incorrect argument supplied:\n"+ args[a] +"\nexiting");
					System.exit(1);
			}
		}
	}
}
