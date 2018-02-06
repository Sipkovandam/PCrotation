package MatrixScripts;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

import RowAnalyses.RowAverageCalculator;
import RowAnalyses.RowCenterer;
import RowAnalyses.RowJobExecutor;
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
	int nThreads=1;
	
	public void run()
	{
		try
		{
			init();
				
			RowJobExecutor rowJobExecutor = new RowJobExecutor(nThreads);
			rowJobExecutor.setWriteFolder(new File(this.writeFn).getParent());
			
			//rowJobs
			RowAverageCalculator rowAverageCalculator = new RowAverageCalculator();
			rowAverageCalculator.setWriteFn(new File(this.averagesWriteFn).getName());
			
			RowCenterer center = new RowCenterer();
			center.setWriteFn(new File(this.writeFn).getName());

			rowJobExecutor.addJob(rowAverageCalculator);
			rowJobExecutor.addJob(center);

			RowJobExecutor.useExecutorsOnFile(rowJobExecutor, fileName);
			
			System.out.println("Done! File written to: " + writeFn);
		}catch(Exception e){e.printStackTrace();}
	}
	
	private void init()
	{
		if(writeFn==null)
		{
			writeFn=FileUtils.addBeforeExtention(this.fileName,"_rowCentered");
		}
		if(averagesWriteFn==null)
		{
			averagesWriteFn=FileUtils.addBeforeExtention(this.fileName,"_averages");
		}
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
