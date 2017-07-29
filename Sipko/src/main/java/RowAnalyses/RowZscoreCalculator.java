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
	
//	final String zScoresName = "zScores";
//
//	RowZscoreCalculator()
//	{
//		this.setValueNames(new String[]{zScoresName});
//	}
	
	@Override
	public void execute(RowJobExecutor rowExecutor, int lineNumber)
	{
		try
		{		
			if(executorWithAvgStdevsToUse == null)
				executorWithAvgStdevsToUse = rowExecutor;
			
			double avg = 0;
			if(rowAverageCalculator!= null && executorWithAvgStdevsToUse.getJobVarIsCalulated(rowAverageCalculator.getAverageName()))
				avg = executorWithAvgStdevsToUse.getJobValue(rowAverageCalculator.getAverageName());
			else
				avg = getAvg(rowExecutor);
			double stdev = 0;
			if(rowStdevCalculator!=null && executorWithAvgStdevsToUse.getJobVarIsCalulated(rowStdevCalculator.getStdevName()))
				stdev = executorWithAvgStdevsToUse.getJobValue(rowStdevCalculator.getStdevName());
			else
				stdev= getStdev(rowExecutor);
			
			double[] zScores = convertToZscores(rowExecutor.getRowName(), rowExecutor.getValues(), avg, stdev);
			writeResult(zScores, rowExecutor, lineNumber);
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

	private void writeResult(double[] zScores, RowJobExecutor rowExecutor, int lineNumber) throws FileNotFoundException, IOException
	{
		String writeLine = rowExecutor.getRowName();
		writeLine=writeLine.concat(FileUtils.doubleArrayToWriteString(zScores));
		writeLine=writeLine.concat("\n");
		super.writeLine(lineNumber, writeLine, rowExecutor, false);
	}


	private double getStdev(RowJobExecutor rowExecutor)
	{
		double variance = org.apache.commons.math3.stat.StatUtils.variance(rowExecutor.getValues());
		double stdev =  java.lang.Math.pow(variance,0.5);
		return stdev;
	}

	private double getAvg(RowJobExecutor rowExecutor)
	{
		double avg = org.apache.commons.math3.stat.StatUtils.mean(rowExecutor.getValues());
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
