package RowAnalyses;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;

import Tools.FileUtils;

public class RowJobExecutorThread implements Runnable
{
	String line;
	int lineNumber;
	ArrayList<RowJobExecutor> rowJobExecutors = null;
	int threadNumber;
	
	boolean[] availableThreadNumbers;
	
	public RowJobExecutorThread(String line,
								int lineNumber,
								ArrayList<RowJobExecutor> rowJobExecutors, boolean[] availableThreadNumbers, int threadNumber)
	{
		this.line = line;
		this.lineNumber = lineNumber;
		this.rowJobExecutors = rowJobExecutors;
		this.threadNumber = threadNumber;
		this.availableThreadNumbers = availableThreadNumbers;
		this.threadNumber = threadNumber;
	}

	@Override
	public void run()
	{
		try
		{
			double startTime = System.nanoTime();
			String rowName = line.split("\t",
										2)[0];
			String[] valuesString = line.split(	"\t",
												2)[1].split("\t");
			
			startTime=printTime(startTime, "SplitTime = ");

			double[] values = getInputValues(valuesString);
			
			startTime=printTime(startTime, "ParseTime = ");
			setRowJobExecutorVariables(	rowJobExecutors,
										rowName,
										values,
										lineNumber);
			
			startTime=printTime(startTime, "SetVariableTime = ");
	
			executeJobs(rowJobExecutors);
			
			startTime=printTime(startTime, "ExecuteTime = ");
			
			this.availableThreadNumbers[threadNumber]=true;
		}catch(Exception e){e.printStackTrace();}
	}

	private double printTime(double startTime, String printString)
	{
//		double time = (System.nanoTime()-startTime)/1000;
//		if(time>3000000)
//			System.out.println(printString + time + "\t threadNumber = " + this.threadNumber + " \t");
//
//		return System.nanoTime();
		return 0;
	}

	private double[] getInputValues(String[] valuesString)
	{
		double[] values =rowJobExecutors.get(0).threadVars[this.threadNumber].getInputValuesAll();
		if(values == null)
			values = FileUtils.convertToDoubleArray(valuesString);
		else
			FileUtils.convertToDoubleArray(valuesString, values);
		return values;
	}

	private void setRowJobExecutorVariables(ArrayList<RowJobExecutor> rowJobExecutors,
											String rowName,
											double[] values,
											int lineNumber) throws FileNotFoundException, IOException
	{
		for (RowJobExecutor rowJobExecutor : rowJobExecutors)
		{
			rowJobExecutor.setRowName(rowName, this.threadNumber);
			
			rowJobExecutor.initiateStorageVariables(this.threadNumber);
			setValuesUsingIncludeIndexes(values, rowJobExecutor);
			
		}
	}
	public void setValuesUsingIncludeIndexes(double[] values, RowJobExecutor rowJobExecutor) throws FileNotFoundException, IOException
	{
		int[] includeIndexes = rowJobExecutor.getIncludeIndexes();
		double[] includeValues = rowJobExecutor.threadVars[this.threadNumber].getIncludeIndexValues();
		if(includeValues == null)
			includeValues=new double[includeIndexes.length];
		
		for (int i = 0; i < includeIndexes.length; i++)
		{
			int e = includeIndexes[i];
			includeValues[i] = values[e];
		}
		rowJobExecutor.threadVars[this.threadNumber].setIncludeIndexValues(includeValues);
	}

	private void executeJobs(ArrayList<RowJobExecutor> rowJobExecutors)
	{
		for (RowJobExecutor rowJobExecutor : rowJobExecutors)
		{
			rowJobExecutor.execute(this.lineNumber, this.threadNumber);
		}
	}
}
