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
	public void execute(RowJobExecutor rowExecutor, int lineNumber, int threadNumber)
	{
		try
		{
			double[] values = rowExecutor.getInputValues(threadNumber);
			double max =  org.apache.commons.math3.stat.StatUtils.max(values);
			super.writeResult(max, rowExecutor, lineNumber, threadNumber);
		}catch(Exception e){e.printStackTrace();}
	}
}
