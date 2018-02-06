package RowAnalyses;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;

import Tools.FileUtils;

public class RowCenterer extends RowJob
{
	private RowAverageCalculator rowAverageCalculator = null;
	private RowJobExecutor executorWithAvgStdevsToUse =null;

//	RowZscoreCalculator()
//	{
//		this.setValueNames(new String[]{zScoresName});
//	}
	
	public RowCenterer()
	{
		this.setHasSingleColHeader(false);
	}
	
	@Override
	public void execute(RowJobExecutor rowExecutor, int lineNumber, int threadNumber)
	{
		try
		{		
			//lo2
			double[] centered = log2( rowExecutor.getInputValues(threadNumber));
			super.writeResult(centered, rowExecutor, lineNumber, threadNumber);
		}catch(Exception e){e.printStackTrace();}
	}
	private double[] log2(double[] values) throws FileNotFoundException, IOException
	{
		double[] log2 = new double[values.length];
		double logVal = Math.log(2);
		for(int v = 0; v < values.length; v++)
		{
			log2[v]= Math.log(values[v] + 1) / logVal;		
		}
		return log2;
	}

	private double getAvg(RowJobExecutor rowExecutor, int threadNumber)
	{
		double avg = org.apache.commons.math3.stat.StatUtils.mean(rowExecutor.getInputValues(threadNumber));
		return avg;
	}
	public RowAverageCalculator getRowAverageCalculator()
	{
		return rowAverageCalculator;
	}
	public void setRowAverageCalculator(RowAverageCalculator rowAverageCalculator)
	{
		this.rowAverageCalculator = rowAverageCalculator;
	}
	public RowJobExecutor getExecutorWithAvgStdevsToUse()
	{
		return executorWithAvgStdevsToUse;
	}
	public void setExecutorWithAvgStdevsToUse(RowJobExecutor executorWithAvgStdevsToUse)
	{
		this.executorWithAvgStdevsToUse = executorWithAvgStdevsToUse;
	}

}
