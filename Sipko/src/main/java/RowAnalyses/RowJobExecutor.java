package RowAnalyses;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import Tools.FileUtils;

public class RowJobExecutor
{
	private String writeFolder = null;
	private List<RowJob> rowJobs = null;
	private int[] includeIndexes = null;
	private ArrayList<String> includeColHeaders = null;
	private String[] dataColHeaders = null;
	
	HashMap<String, Integer> valueNameToJobVarIndex = null;

	ThreadVars[] threadVars = null; 
	
//	HashMap<String,ArrayList<String>> resultSlots = new HashMap<String, ArrayList<String>>(); //something for making it multithreaded
//	HashMap<String,Integer> resultsWrittenIndex = new HashMap<String,Integer>();

	//lineWriter for each job
	HashMap<String,JobLineWriter> lineWriters = new HashMap<String,JobLineWriter>();

	
	public RowJobExecutor(int nThreads)
	{
		super();
		threadVars = new ThreadVars[nThreads];
		for(int t =0; t < nThreads; t++)
		{
			threadVars[t] = new ThreadVars();
		}
	}

	public void execute(int lineNumber, int threadNumber)
	{
		for(RowJob rowJob : rowJobs)
		{
			rowJob.execute(this, lineNumber, threadNumber);
		}
	}

	public void addJob(RowJob rowJob)
	{
		if(rowJobs == null)
			rowJobs = new ArrayList<RowJob>();
		rowJobs.add(rowJob);
	}

	public void setJobs(List<RowJob> rowJobs)
	{
		this.rowJobs=rowJobs;
	}

	public void executeHeaderJob(String header)
	{
		this.dataColHeaders=executeHeaderJobForWriteFolder(header);
	}
	
	public String[] executeHeaderJobForWriteFolder(String header)
	{
		String[] dataHeaders = header.split("\t",2)[1].split("\t");
		
		if(this.includeIndexes ==null)
			setIncludeIndexes(dataHeaders);
			
		String[] dataColHeaders = new String[this.includeIndexes.length];
		for(int i = 0; i < this.includeIndexes.length; i++)
		{
			int e = this.includeIndexes[i];
			dataColHeaders[i]=dataHeaders[e];
		}		
		return dataColHeaders;
	}
//
//	public void setValuesUsingIncludeIndexes(double[] values, int threadNumber)
//	{
//		ThreadVars threadVars = this.threadVars[threadNumber];
//		double[] tempValues = threadVars.getIncludeIndexValues();
//		
//		if(tempValues==null)
//			tempValues=new double[includeIndexes.length];
//		
//		for(int i = 0; i < includeIndexes.length; i++)
//		{
//			int e = includeIndexes[i];
//			tempValues[i]=values[e];
//
//			for(RowJob rowJob : rowJobs)
//			{
//				rowJob.executeOnInitiation(this, tempValues[i], threadNumber);
//			}
//		}
//		threadVars.setIncludeIndexValues(tempValues);
//	}

	public int[] getIncludeIndexes()
	{
		return includeIndexes;
	}

	public void setIncludeIndexes(int[] includeIndexes)
	{
		this.includeIndexes = includeIndexes;
	}

	public String[] getDataColHeaders()
	{
		return dataColHeaders;
	}

	public void setIncludeIndexes(String[] dataHeaders)
	{
		if(this.getIncludeIndexes() == null && this.includeColHeaders != null)
		{		
			HashMap<String, Integer> colHeaderToColIndex = FileUtils.arrayToHashMap(dataHeaders);
			
			int nPresent = 0;
			for(String includeHeader : this.includeColHeaders)
			{
				if(colHeaderToColIndex.get(includeHeader) != null)
					nPresent++;
				else
				{
//					for(String colHeader: colHeaderToColIndex.keySet())
//						System.out.println("colheader = \t" + colHeader);
					System.out.println("RowJobExecutorWarning - IncludeHeader:" + includeHeader + ", is missing in the columnHeaders and is not included");
				}
			}
			this.includeIndexes = new int[nPresent];
			int n =0;
			for(String includeHeader : this.includeColHeaders)
			{
				if(colHeaderToColIndex.containsKey(includeHeader))
				{
					int col = colHeaderToColIndex.get(includeHeader);
					this.includeIndexes[n] = col;
					n++;
				}
			}
		}
		else
		{
			System.out.println("DataColumns to include not initiated");
		}
	}

	public String getWriteFolder()
	{
		return writeFolder;
	}

	public void setWriteFolder(String writeFolder)
	{
		this.writeFolder = writeFolder;
	}

	public void setIncludeColHeaders(ArrayList<String> includeColnames)
	{
		this.includeColHeaders = includeColnames;
	}

	public List<RowJob> getRowJobs()
	{
		return rowJobs;
	}

	public void setRowJobs(List<RowJob> rowJobs)
	{
		this.rowJobs = rowJobs;
	}
	
	public void initiateStorageVariables(int threadNumber)
	{
		if(valueNameToJobVarIndex==null)
		{
			valueNameToJobVarIndex= new HashMap<String, Integer>();
			
			for(RowJob rowJob : rowJobs)
			{
				String[] valueNames = rowJob.getValueNames();
				
				if(valueNames!= null)
					for(String valueName:valueNames)
						valueNameToJobVarIndex.put(valueName, valueNameToJobVarIndex.size());
			}
		}

		if(this.threadVars[threadNumber].getThread_JobVars()==null)
			this.threadVars[threadNumber].setThread_JobVars(new double[valueNameToJobVarIndex.size()]);
		if(this.threadVars[threadNumber].getThread_JobVarIsCalulated()==null)
			this.threadVars[threadNumber].setThread_JobVarIsCalulated(new boolean[valueNameToJobVarIndex.size()]);
		else
			this.threadVars[threadNumber].setThread_JobVarIsCalulatedFalse();
	}

	public void writeHeaders() throws FileNotFoundException, IOException
	{
		for(RowJob rowJob : rowJobs)
		{
			JobLineWriter jobLineWriter= new JobLineWriter(this.writeFolder, rowJob.getWriteFn(), this.getDataColHeaders(), rowJob.hasSingleColHeader());
			lineWriters.put(rowJob.getWriteFn(), jobLineWriter);
		}
		
	}

	public double getJobValue(String valueName, int threadNumber)
	{
		int valueNameIndex = this.valueNameToJobVarIndex.get(valueName);
		return this.threadVars[threadNumber].getThread_JobVars()[valueNameIndex];
	}

	public void setJobValue(String valueName, double value, int threadNumber)
	{
		int valueIndex = this.valueNameToJobVarIndex.get(valueName);	
		this.threadVars[threadNumber].setThread_JobVars(valueIndex, value);
		this.threadVars[threadNumber].setThread_JobVarIsCalulated(valueIndex,true);
	}

	public HashMap<String, Integer> getValueNameToJobVarIndex()
	{
		return valueNameToJobVarIndex;
	}

	public void setValueNameToJobVarIndex(HashMap<String, Integer> valueNameToJobVarIndex)
	{
		this.valueNameToJobVarIndex = valueNameToJobVarIndex;
	}

	public boolean getJobVarIsCalulated(String valueName, int threadNumber)
	{
		int valueIndex = this.valueNameToJobVarIndex.get(valueName);	
		return this.threadVars[threadNumber].getThread_JobVarIsCalulated(valueIndex);
	}

	public void write(String writeLine, String writeFn, boolean hasSingleColHeader, int lineNumber,int threadNumber) throws FileNotFoundException, IOException, InterruptedException
	{
		JobLineWriter writer = getJobLineWriter(writeFn,hasSingleColHeader);
		writer.writeIntoBufferAndWriteBuffer(lineNumber,writeLine, threadNumber);
	}

	private JobLineWriter getJobLineWriter(	String writeFn,
											boolean hasSingleColHeader) throws FileNotFoundException, IOException
	{
		JobLineWriter jobLineWriter = lineWriters.get(writeFn);
		return jobLineWriter;
	}


	public void closeWriters() throws IOException
	{
		for(JobLineWriter jobLineWriter: lineWriters.values())
		{
			jobLineWriter.close();
		}
	}

	public void setRowName(	String rowName,
							int threadNumber)
	{
		this.threadVars[threadNumber].setRowName(rowName);
	}

	public double[] getInputValues(int threadNumber)
	{
		return this.threadVars[threadNumber].getIncludeIndexValues();
	}

	public String getRowName(int threadNumber)
	{
		return this.threadVars[threadNumber].getRowName();
	}

	public HashMap<String, JobLineWriter> getLineWriters()
	{
		return lineWriters;
	}
}
