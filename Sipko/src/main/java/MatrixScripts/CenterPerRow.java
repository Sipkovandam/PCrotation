package MatrixScripts;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

import Tools.FileUtils;
import Tools.Script;

public class CenterPerRow extends Script<CenterPerRow>
{
	String fileNameComment = "/root/directory/valueMatrix.txt; OPTIONAL; Input file for which the rowAverages should be calculated";
	String fileName = null;
	String averagesWriteFnComment = "/root/directory/rowAverages.txt; OPTIONAL; File where the averages per row should be written";
	String averagesWriteFn = null;
	String writeFnComment = "/root/directory/averges.txt; OPTIONAL; Path where the file containing the row averages and standard deviations will be written";
	String writeFn = null;
	
	public void run()
	{
		try
		{
			if(writeFn == null)
				writeFn = FileUtils.removeExtention(fileName)+"_centered.txt";
			if(averagesWriteFn == null)
				averagesWriteFn = FileUtils.removeExtention(fileName)+"_rowAverages.txt";
			
			MyMatrix matrix = new MyMatrix(fileName);
			MyMatrix rowAverages = matrix.calcAvgRows();
			System.out.println(rowAverages.values[0][0]);
			matrix.adjustForAverageAllrows(rowAverages);

			rowAverages.write(averagesWriteFn);
			matrix.write(writeFn);
			
			System.out.println("Done! File written to: " + writeFn);
		}catch(Exception e){e.printStackTrace();}
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

	public String getAveragesWriteFn()
	{
		return averagesWriteFn;
	}

	public void setAveragesWriteFn(String averagesWriteFn)
	{
		this.averagesWriteFn = averagesWriteFn;
	}
}
