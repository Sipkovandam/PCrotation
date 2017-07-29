package RowAnalyses;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.apache.commons.math3.stat.descriptive.rank.Median;

import Tools.FileUtils;

public class RowColumnGetter extends RowJob
{
	@Override
	public void execute(RowJobExecutor rowExecutor, int lineNumber)
	{
		try
		{
			double[] values = rowExecutor.getValues();
			writeResult(values, rowExecutor, lineNumber);
		}catch(Exception e){e.printStackTrace();}
	}

	private void writeResult(double[] values, RowJobExecutor rowExecutor, int lineNumber) throws FileNotFoundException, IOException
	{
		String writeLine = rowExecutor.getRowName().concat(FileUtils.doubleArrayToWriteString(values)).concat("\n");
		super.writeLine(lineNumber, writeLine, rowExecutor, false);
	}
}
