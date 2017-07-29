package RowAnalyses;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.apache.commons.math3.stat.descriptive.rank.Median;

import Tools.FileUtils;

public class RowMaxGetter extends RowJob
{
	final String maxName = "max";
	
	public RowMaxGetter()
	{
		this.setValueNames(new String[]{maxName});
	}
	
	@Override
	public void execute(RowJobExecutor rowExecutor, int lineNumber)
	{
		try
		{
			double max =rowExecutor.getJobValue(maxName);
			writeResult(max, rowExecutor, lineNumber);
		}catch(Exception e){e.printStackTrace();}
	}
	
	@Override
	public void executeOnInitiation(RowJobExecutor rowExecutor, double value)
	{
		if(rowExecutor.getJobValue(maxName)<value)
			rowExecutor.setJobValue(maxName,value);
	}

	private void writeResult(double max, RowJobExecutor rowExecutor, int lineNumber) throws FileNotFoundException, IOException
	{
		String writeLine = rowExecutor.getRowName().concat("\t").concat(Double.toString(max).concat("\n"));
		super.writeLine(lineNumber, writeLine, rowExecutor, true);
	}
	
//	public double getMax(double[] values)
//	{
//		double max = Double.NEGATIVE_INFINITY;
//		for(double value:values)
//		{
//			if(value>max)
//				max=value;
//		}
//		return max;
//	}
	
	
}
