package RowAnalyses;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.apache.commons.math3.stat.descriptive.rank.Median;

import Tools.FileUtils;

public class RowAboveCutoffCounter extends RowJob
{
	double cutoff = 0;
	boolean absolute = false;
	
	final String aboveCutoffCountName = "aboveCutoffCount";
	
	public RowAboveCutoffCounter()
	{
		this.setValueNames(new String[]{aboveCutoffCountName});
	}
	
	@Override
	public void execute(RowJobExecutor rowExecutor, int lineNumber)
	{
		try
		{
			double aboveCutoffCount = rowExecutor.getJobValue(aboveCutoffCountName);
			
			writeResult(aboveCutoffCount, rowExecutor, lineNumber);
		}catch(Exception e){e.printStackTrace();}
	}
	
	@Override
	public void executeOnInitiation(RowJobExecutor rowExecutor, double value)
	{
		if(absolute)
		{
			if(Math.abs(rowExecutor.getJobValue(aboveCutoffCountName))>cutoff)
				rowExecutor.setJobValue(aboveCutoffCountName,rowExecutor.getJobValue(aboveCutoffCountName)+1);
		}
		else
			if(rowExecutor.getJobValue(aboveCutoffCountName)>cutoff)
				rowExecutor.setJobValue(aboveCutoffCountName,rowExecutor.getJobValue(aboveCutoffCountName)+1);
	}

	private void writeResult(double aboveCutoffCount, RowJobExecutor rowExecutor, int lineNumber) throws FileNotFoundException, IOException
	{
		String writeLine = rowExecutor.getRowName().concat("\t").concat(Integer.toString(((int)aboveCutoffCount)).concat("\n"));
		super.writeLine(lineNumber, writeLine, rowExecutor, true);
	}

	public double getCutoff()
	{
		return cutoff;
	}

	public void setCutoff(double cutoff)
	{
		this.cutoff = cutoff;
	}

	public boolean isAbsolute()
	{
		return absolute;
	}

	public void setAbsolute(boolean absolute)
	{
		this.absolute = absolute;
	}

}
