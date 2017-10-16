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
	public void execute(RowJobExecutor rowExecutor, int lineNumber, int threadNumber)
	{
		try
		{
			double[] values = rowExecutor.getInputValues(threadNumber);
			int aboveCutoffCount = getAboveCutoffCount(values, cutoff, absolute);
			
			rowExecutor.setJobValue(aboveCutoffCountName,aboveCutoffCount, threadNumber);
			
			super.writeResult(aboveCutoffCount, rowExecutor, lineNumber, threadNumber);
		}catch(Exception e){e.printStackTrace();}
	}
	
	private int getAboveCutoffCount(	double[] values,
										double cutoff, boolean absolute)
	{
		int count = 0;
		for(double value:values)
		{
			if(absolute)
				value = Math.abs(value);
			if(value>cutoff)
				count++;
		}
		return count;
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
