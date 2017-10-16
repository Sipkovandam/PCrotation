package RowAnalyses;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;

import Tools.FileUtils;

public class RowAverageCalculator extends RowJob
{
	final String averageName = "average";
	public RowAverageCalculator()
	{
		this.setValueNames(new String[]{averageName});
	}
	
	@Override
	public void execute(RowJobExecutor rowExecutor, int lineNumber, int threadNumber)
	{
		try
		{
			double[] values = rowExecutor.getInputValues(threadNumber);
			double avg = org.apache.commons.math3.stat.StatUtils.mean(values);
			rowExecutor.setJobValue(averageName, avg, threadNumber);
			super.writeResult(avg, rowExecutor, lineNumber, threadNumber);
		}catch(Exception e){e.printStackTrace();}
	}

	public String getAverageName()
	{
		return averageName;
	}
}

