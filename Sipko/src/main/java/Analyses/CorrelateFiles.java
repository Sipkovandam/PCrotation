package Analyses;

import Tools.FileUtils;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

import MatrixScripts.MyMatrix;

public class CorrelateFiles 
{
	static String fn1 = "E:/Groningen/Data/Juha/Genes31995/31.07.pc1.illumina.genes.expressed.QuantileNorm/31.07.pc1.illumina.genes.expressed.QuantileNormalized.Log2Transformed_SMALL.txt";
	static String fn2 = "E:/Groningen/Data/Juha/Genes31995/31.07.pc1.illumina.genes.expressed.QuantileNorm/31.07.pc1.illumina.genes.expressed.QuantileNormalized.Log2Transformed_SMALL_test.txt";
	static String writeName1 = null;
	static String writeName2 = null;

	public static void main (String[] args)
	{
		checkArgs(args);
		if(writeName1 == null)
			writeName1 = FileUtils.replaceEnd(fn1, "_VSfile2CorrRows.txt");
		if(writeName2 == null)
			writeName2 = FileUtils.replaceEnd(fn1, "_VSfile2CorrColumns.txt");
		
		
		MyMatrix m1 = new MyMatrix(fn1);
		MyMatrix m2 = new MyMatrix(fn2);

		m1.keepRows(m2);//keep only rows present in both files and put them in the same order
		
		MyMatrix corrRows = getCorPerRow(m1, m2);
		corrRows.write(writeName1);
		System.out.println("File written to:" + writeName1);
		m1.transpose();
		m2.transpose();
		m1.keepRows(m2);
		
		MyMatrix corrCols = getCorPerRow(m1, m2);

		corrCols.write(writeName2);
		System.out.println("File written to:" + writeName2);
	}
	private static MyMatrix getCorPerRow(MyMatrix m1, MyMatrix  m2) {
		
		PearsonsCorrelation correlator = new PearsonsCorrelation();
		MyMatrix corrRows = new MyMatrix(m1.rows(), 1);
		corrRows.rowNames= m1.rowNames;
		corrRows.colNames= new String[]{"Correlation"};
		for(int r = 0; r < m1.rows(); r++)
		{
			double cor = correlator.correlation(m1.values[r], m2.values[r]);
			corrRows.matrix.set(r, 0, cor);
		}
		return corrRows;
	}
	static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko") && args.length < 1)
			return;
		if(args.length < 1)
		{
			System.out.println("Script requires the following argumetns:\n"
					+ "1. fn1=<filename1.txt> - First file to compare\n"
					+ "2. fn2=<filename2.txt> - Second file to compare (to the first file)\n"
					+ "3. writeFN1=<writeFilename1.txt> - Name of the row correlation file to be written (default=<FN1>+_VSfile2CorrRows.txt\n"
					+ "4. writeFN2=<writeFilename2.txt> - Name of the column correlation file to be written (default=<FN1>+_VSfile2CorrColumns.txt");
			System.exit(1);
		}
		
		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase())
			{
//				var = new JSONutil<Vars>().read(var.JSON_FN, var);
				case "fn1":
					fn1 =value;
					break;
				case "fn2":
					fn2 =value;
					break;
				case "writefn1":
					writeName1 =value;
					break;
				case "writefn2":
					writeName2 =value;
					break;
				default:
					System.out.println("Incorrect argument supplied:\n"+ args[a] +"\nexiting");
					System.exit(1);
			}
		}
	}
}
