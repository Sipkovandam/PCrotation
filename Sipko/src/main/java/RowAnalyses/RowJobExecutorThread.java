//package RowAnalyses;
//
//import java.io.BufferedWriter;
//import java.util.ArrayList;
//import java.util.HashMap;
//
//import Tools.FileUtils;
//
//public class RowJobExecutorThread implements Runnable
//{
//	String line;
//	int lineNumber;
//	ArrayList<RowJobExecutor> rowJobExecutors = null;
//	
//	public RowJobExecutorThread(String line,
//								int lineNumber,
//								ArrayList<RowJobExecutor> rowJobExecutors)
//	{
//		this.line = line;
//		this.lineNumber = lineNumber;
//		this.rowJobExecutors = rowJobExecutors;
//	}
//
//	@Override
//	public void run()
//	{
//		String rowName = line.split("\t",
//									2)[0];
//		String[] valuesString = line.split(	"\t",
//											2)[1].split("\t");
//		double[] values = FileUtils.convertToDoubleArray(valuesString);
//		
//		setRowJobExecutorVariables(	rowJobExecutors,
//									rowName,
//									values,
//									lineNumber);
//
//		executeJobs(rowJobExecutors);
//	}
//
//	private void setRowJobExecutorVariables(ArrayList<RowJobExecutor> rowJobExecutors,
//											String rowName,
//											double[] values,
//											int lineNumber)
//	{
//		for (RowJobExecutor rowJobExecutor : rowJobExecutors)
//		{
//			rowJobExecutor.setRowName(rowName);
//			
//			setValuesUsingIncludeIndexes(values, rowJobExecutor);
//			
//			for(RowJob rowJob : rowJobExecutor.getRowJobs())
//			{
//				//rowJob.addResultSlots(lineNumber);
//			}
//				
//		}
//	}
//	public void setValuesUsingIncludeIndexes(double[] values, RowJobExecutor rowJobExecutor)
//	{
//		int[] includeIndexes = rowJobExecutor.getIncludeIndexes();
//		double[] includeValues= new double[includeIndexes.length];
//		for(int i = 0; i < includeIndexes.length; i++)
//		{
//			int e = includeIndexes[i];
//			includeValues[i]=values[e];
//		}
//		rowJobExecutor.setValues(includeValues);
//	}
//
//	private void executeJobs(ArrayList<RowJobExecutor> rowJobExecutors)
//	{
//		for (RowJobExecutor rowJobExecutor : rowJobExecutors)
//		{
//			rowJobExecutor.execute(this.lineNumber);
//		}
//	}
//}
