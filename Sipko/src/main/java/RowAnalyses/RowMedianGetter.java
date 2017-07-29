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
	public void execute(RowJobExecutor rowExecutor, int lineNumber)
	{
		try
		{
			double median =getMedian(rowExecutor.getValues());
			rowExecutor.setJobValue(medianName, median);
			writeResult(median, rowExecutor, lineNumber);
		}catch(Exception e){e.printStackTrace();}
	}

	private void writeResult(double avg, RowJobExecutor rowExecutor, int lineNumber) throws FileNotFoundException, IOException
	{
		String writeLine = rowExecutor.getRowName().concat("\t").concat(Double.toString(avg).concat("\n"));
		super.writeLine(lineNumber, writeLine, rowExecutor, true);
	}
	
	public double getMedian(double[] values)
	{
		Median median = new Median();
		double medianValue = median.evaluate(values);
		return medianValue;
	}
}
