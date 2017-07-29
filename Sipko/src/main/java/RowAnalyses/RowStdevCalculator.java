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
	public void execute(RowJobExecutor rowExecutor, int lineNumber)
	{
		try
		{
			double variance = org.apache.commons.math3.stat.StatUtils.variance(rowExecutor.getValues());
			double stdev =  java.lang.Math.pow(variance,0.5);
			
			rowExecutor.setJobValue(stdevName, stdev);
			
			writeResult(stdev, rowExecutor, lineNumber);
		}catch(Exception e){e.printStackTrace();}
	}

	private void writeResult(double avg, RowJobExecutor rowExecutor, int lineNumber) throws FileNotFoundException, IOException
	{
		String writeLine = rowExecutor.getRowName().concat("\t").concat(Double.toString(avg).concat("\n"));
		super.writeLine(lineNumber, writeLine, rowExecutor, true);
	}
	public String getStdevName()
	{
		return stdevName;
	}
}
