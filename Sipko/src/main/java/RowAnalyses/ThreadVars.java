package RowAnalyses;

import java.util.Arrays;

public class ThreadVars
{
	private String rowName = null;
	private double[] inputValuesAll = null;//shared across all threads and executors (for 1 file)
	private double[] includeIndexValues = null;
	private double[] thread_JobVars = null;
	private boolean[] thread_JobVarIsCalulated = null;
	private int threadNumber = 0;
	
	public String getRowName()
	{
		return rowName;
	}
	public void setRowName(String rowName)
	{
		this.rowName = rowName;
	}
	public double[] getIncludeIndexValues()
	{
		return includeIndexValues;
	}
	public void setIncludeIndexValues(double[] inputValues)
	{
		this.includeIndexValues = inputValues;
	}
	public double[] getThread_JobVars()
	{
		return thread_JobVars;
	}
	public void setThread_JobVars(double[] thread_JobVars)
	{
		this.thread_JobVars = thread_JobVars;
	}
	public boolean[] getThread_JobVarIsCalulated()
	{
		return thread_JobVarIsCalulated;
	}
	public void setThread_JobVarIsCalulated(boolean[] thread_JobVarIsCalulated)
	{
		this.thread_JobVarIsCalulated = thread_JobVarIsCalulated;
	}
	public int getThreadNumber()
	{
		return threadNumber;
	}
	public void setThreadNumber(int threadNumber)
	{
		this.threadNumber = threadNumber;
	}
	public void setThread_JobVars(	int valueIndex,
									double value)
	{
		this.thread_JobVars[valueIndex] = value;
	}
	public void setThread_JobVarIsCalulated(int valueIndex,
											boolean value)
	{
		this.thread_JobVarIsCalulated[valueIndex] = value;
	}
	public boolean getThread_JobVarIsCalulated(int valueIndex)
	{
		return this.thread_JobVarIsCalulated[valueIndex] ;
	}
	public double[] getInputValuesAll()
	{
		return inputValuesAll;
	}
	public void setInputValuesAll(double[] inputValuesAll)
	{
		this.inputValuesAll = inputValuesAll;
	}
	public void setThread_JobVarIsCalulatedFalse()
	{
		Arrays.fill(this.thread_JobVarIsCalulated, false);	
	}
	
}
