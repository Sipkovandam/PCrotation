package RowAnalyses;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.apache.commons.math3.stat.descriptive.rank.Median;

import Tools.FileUtils;

public class RowColumnGetter extends RowJob
{
	public RowColumnGetter()
	{
		this.setHasSingleColHeader(false);
	}
	
	@Override
	public void execute(RowJobExecutor rowExecutor, int lineNumber, int threadNumber)
	{
		try
		{
			double[] values = rowExecutor.getInputValues(threadNumber);
			super.writeResult(values, rowExecutor, lineNumber, threadNumber);
		}catch(Exception e){e.printStackTrace();}
	}
}
