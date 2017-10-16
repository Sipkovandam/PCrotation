package RowAnalyses;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.apache.commons.math3.stat.descriptive.rank.Median;

public class RowStdevCalculator extends RowJob
{
	final String stdevName = "stdev";
	
	public RowStdevCalculator()
	{
		this.setValueNames(new String[]{stdevName});
	}
	
	@Override
	public void execute(RowJobExecutor rowExecutor, int lineNumber, int threadNumber)
	{
		try
		{
			double variance = org.apache.commons.math3.stat.StatUtils.variance(rowExecutor.getInputValues(threadNumber));
			double stdev =  java.lang.Math.pow(variance,0.5);
			
			rowExecutor.setJobValue(stdevName, stdev, threadNumber);
			
			super.writeResult(stdev, rowExecutor, lineNumber, threadNumber);
		}catch(Exception e){e.printStackTrace();}
	}

	public String getStdevName()
	{
		return stdevName;
	}
}
