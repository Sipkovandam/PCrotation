package MatrixScripts;

import Tools.FileUtils;
import Tools.Script;

public class Transpose extends Script<Transpose>
{
	/**
	 * 
	 */
	private static final long serialVersionUID = -3861427697209249113L;
	//Transposes the input file.
	String fileName = null;
	String writeFn = null;
	
	public void run()
	{
		//String fileName = "E:/Groningen/Data/PublicSamples/05-2016/22214Samples_Voom_Correl/SAMPLE_NormalizedLog2_10Rows.txt";//
		MyMatrix matrix = new MyMatrix(fileName);
		//matrix.print(5,5);
		System.out.println("Input rows = " + matrix.rowNames.length + " columns = " + matrix.colNames.length);
		
		if(writeFn ==null)
			writeFn=FileUtils.removeExtention(fileName)+ "_transposed.txt.gz";
		matrix.writeTransposed(writeFn);
		
//		matrix.write(fileName.replace(".txt", "_first10.txt"), -1,10);
//		matrix.write(fileName.replace(".txt", "_first100geneEvs.txt"), -1,100);
//		matrix.write(fileName.replace(".txt", "_transposed_100.txt"), 100, -1);
//		matrix.write(fileName.replace(".txt", "_transposed_300.txt"), 300, -1);
//		matrix.write(fileName.replace(".txt", "_transposed_500.txt"), 500, -1);
//		matrix.write(fileName.replace(".txt", "_transposed_1000.txt"), 1000, -1);
		System.out.println("Done, File written at:" + writeFn);
	}

	public String getFileName()
	{
		return fileName;
	}

	public void setFileName(String fileName)
	{
		this.fileName = fileName;
	}

	public String getWriteFn()
	{
		return writeFn;
	}

	public void setWriteFn(String writeFn)
	{
		this.writeFn = writeFn;
	}
}
