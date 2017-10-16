package RowAnalyses;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.apache.commons.math3.stat.descriptive.rank.Median;

import Tools.FileUtils;

public class RowBelowCutoffCounter extends RowJob
{
	double cutoff = 0;
	boolean absolute = false;
	
	final String belowCutoffCountName = "belowCutoffCount";
	
	public RowBelowCutoffCounter()
	{
		this.setValueNames(new String[]{belowCutoffCountName});
	}
	
	@Override
	public void execute(RowJobExecutor rowExecutor, int lineNumber, int threadNumber)
	{
		try
		{
			double[] values = rowExecutor.getInputValues(threadNumber);
			int belowCutoffCount = getBelowCutoffCount(values, cutoff, absolute);
			
			rowExecutor.setJobValue(belowCutoffCountName,belowCutoffCount, threadNumber);
			
			super.writeResult(belowCutoffCount, rowExecutor, lineNumber, threadNumber);
		}catch(Exception e){e.printStackTrace();}
	}
	
	private int getBelowCutoffCount(	double[] values,
										double cutoff, boolean absolute)
	{
		int count = 0;
		for(double value:values)
		{
			if(absolute)
				value = Math.abs(value);
			if(value<cutoff)
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
