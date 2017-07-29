package RowAnalyses;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;

import Tools.FileUtils;

public class RowAverageCalculator extends RowJob
{
	final String averageName = "average";
	final String sumName = "sum";
	
	public RowAverageCalculator()
	{
		this.setValueNames(new String[]{averageName, sumName});
	}
	
	@Override
	public void execute(RowJobExecutor rowExecutor, int lineNumber)
	{
		try
		{
			double avg = rowExecutor.getJobValue(sumName)/rowExecutor.getValues().length;

			rowExecutor.setJobValue(averageName, avg);
			writeResult(avg, rowExecutor, lineNumber);
		}catch(Exception e){e.printStackTrace();}
	}

	@Override
	public void executeOnInitiation(RowJobExecutor rowExecutor, double value)
	{
		rowExecutor.setJobValue(sumName, rowExecutor.getJobValue(sumName)+value);
	}
	
	private void writeResult(double avg, RowJobExecutor rowExecutor, int lineNumber) throws FileNotFoundException, IOException
	{
		String writeLine = rowExecutor.getRowName().concat("\t").concat(Double.toString(avg).concat("\n"));
		super.writeLine(lineNumber, writeLine, rowExecutor, true);
	}

	public String getAverageName()
	{
		return averageName;
	}

	public String getSumName()
	{
		return sumName;
	}
}

