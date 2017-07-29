package RowAnalyses;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import Tools.FileUtils;

public class RowJobExecutor extends Row
{
	private String writeFolder = null;
	private List<RowJob> rowJobs = null;
	private int[] includeIndexes = null;
	private ArrayList<String> includeColHeaders = null;
	private String[] dataColHeaders = null;
	
	HashMap<String, Integer> valueNameToJobVarIndex = null;
	private double[] jobVars = null;
	private boolean[] jobVarIsCalulated = null;
	
//	HashMap<String,ArrayList<String>> resultSlots = new HashMap<String, ArrayList<String>>(); //something for making it multithreaded
//	HashMap<String,Integer> resultsWrittenIndex = new HashMap<String,Integer>();

	HashMap<String,BufferedWriter> fileWriters = new HashMap<String,BufferedWriter>();
	
	//public HashMap<String,BufferedWriter> fileWriters = new HashMap<String,BufferedWriter>();//a writer for each job
	
	
	public RowJobExecutor()
	{
		super();
	}

	public void execute(int lineNumber)
	{
		for(RowJob rowJob : rowJobs)
		{
			rowJob.execute(this, lineNumber);
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


	

	public void closeWriters() throws IOException
	{
		for(BufferedWriter writer: fileWriters.values())
		{
			writer.close();
		}
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

	public void setValuesUsingIncludeIndexes(double[] values)
	{
		this.values= new double[includeIndexes.length];
		for(int i = 0; i < includeIndexes.length; i++)
		{
			int e = includeIndexes[i];
			this.values[i]=values[e];

			for(RowJob rowJob : rowJobs)
			{
				rowJob.executeOnInitiation(this, this.values[i]);
			}
		}
	}

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

	public HashMap<String, BufferedWriter> getFileWriters()
	{
		return fileWriters;
	}

	public void setFileWriters(HashMap<String, BufferedWriter> fileWriters)
	{
		this.fileWriters = fileWriters;
	}
	
	public void initiateStorageVariables()
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

		this.jobVars = new double[valueNameToJobVarIndex.size()];
		this.jobVarIsCalulated = new boolean[valueNameToJobVarIndex.size()];
	}


	public double getJobValue(String valueName)
	{
		int valueIndex = this.valueNameToJobVarIndex.get(valueName);
		return this.jobVars[valueIndex];
	}

	public void setJobValue(String valueName, double value)
	{
		int valueIndex = this.valueNameToJobVarIndex.get(valueName);	
		this.jobVars[valueIndex]= value;
		this.jobVarIsCalulated[valueIndex]=true;
	}

	public HashMap<String, Integer> getValueNameToJobVarIndex()
	{
		return valueNameToJobVarIndex;
	}

	public void setValueNameToJobVarIndex(HashMap<String, Integer> valueNameToJobVarIndex)
	{
		this.valueNameToJobVarIndex = valueNameToJobVarIndex;
	}

	public boolean getJobVarIsCalulated(String valueName)
	{
		int valueIndex = this.valueNameToJobVarIndex.get(valueName);	
		return this.jobVarIsCalulated[valueIndex];
	}
}
