package RowAnalyses;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;

import Tools.FileUtils;

public class RowZscoreCalculator extends RowJob
{
	private RowAverageCalculator rowAverageCalculator = null;
	private RowStdevCalculator rowStdevCalculator = null;
	
	private RowJobExecutor executorWithAvgStdevsToUse =null;
	
	final String zScoresName = "zScores";

//	RowZscoreCalculator()
//	{
//		this.setValueNames(new String[]{zScoresName});
//	}
	
	public RowZscoreCalculator()
	{
		this.setHasSingleColHeader(false);
	}
	
	@Override
	public void execute(RowJobExecutor rowExecutor, int lineNumber, int threadNumber)
	{
		try
		{		
			if(executorWithAvgStdevsToUse == null)
				executorWithAvgStdevsToUse = rowExecutor;
			
			//calculate averages
			double avg = 0;
			if(rowAverageCalculator!= null && executorWithAvgStdevsToUse.getJobVarIsCalulated(rowAverageCalculator.getAverageName(), threadNumber))
				avg = executorWithAvgStdevsToUse.getJobValue(rowAverageCalculator.getAverageName(), threadNumber);
			else
				avg = getAvg(rowExecutor, threadNumber);
			double stdev = 0;
			//calculate standard deviations
			if(rowStdevCalculator!=null && executorWithAvgStdevsToUse.getJobVarIsCalulated(rowStdevCalculator.getStdevName(), threadNumber))
				stdev = executorWithAvgStdevsToUse.getJobValue(rowStdevCalculator.getStdevName(), threadNumber);
			else
				stdev= getStdev(rowExecutor, threadNumber);
			
			//z-scores
			double[] zScores = convertToZscores(rowExecutor.getRowName(threadNumber), rowExecutor.getInputValues(threadNumber), avg, stdev);
			super.writeResult(zScores, rowExecutor, lineNumber, threadNumber);
		}catch(Exception e){e.printStackTrace();}
	}
	private double[] convertToZscores(String rowName, double[] values, double avg, double stdev) throws FileNotFoundException, IOException
	{
		double[] zScores = new double[values.length];
		for(int v = 0; v < values.length; v++)
		{
			zScores[v]=(values[v]-avg)/stdev;
		}
		return zScores;
	}

	private double getStdev(RowJobExecutor rowExecutor, int threadNumber)
	{
		double variance = org.apache.commons.math3.stat.StatUtils.variance(rowExecutor.getInputValues(threadNumber));
		double stdev =  java.lang.Math.pow(variance,0.5);
		return stdev;
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
	public RowStdevCalculator getRowStdevCalculator()
	{
		return rowStdevCalculator;
	}
	public void setRowStdevCalculator(RowStdevCalculator rowStdevCalculator)
	{
		this.rowStdevCalculator = rowStdevCalculator;
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
