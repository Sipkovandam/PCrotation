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
		MatrixString matrix = new MatrixString(fileName);
		//matrix.print(5,5);
		System.out.println("Input rows = " + matrix.rowNames.length + " columns = " + matrix.colNames.length);
		
		if(writeFn ==null)
			writeFn=FileUtils.removeExtention(fileName)+ "_transposed.txt.gz";
		matrix.writeTransposed(writeFn);
		
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
