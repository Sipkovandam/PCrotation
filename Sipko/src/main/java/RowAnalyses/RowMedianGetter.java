package RowAnalyses;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.apache.commons.math3.stat.descriptive.rank.Median;

import Tools.FileUtils;

public class RowMedianGetter extends RowJob
{
	final String medianName = "median";
	
	public RowMedianGetter()
	{
		this.setValueNames(new String[]{medianName});
	}
	
	@Override
	public void execute(RowJobExecutor rowExecutor, int lineNumber, int threadNumber)
	{
		try
		{
			double median =getMedian(rowExecutor.getInputValues(threadNumber));
			rowExecutor.setJobValue(medianName, median, threadNumber);
			super.writeResult(median, rowExecutor, lineNumber, threadNumber);
		}catch(Exception e){e.printStackTrace();}
	}
	
	public double getMedian(double[] values)
	{
		Median median = new Median();
		double medianValue = median.evaluate(values);
		return medianValue;
	}
}
