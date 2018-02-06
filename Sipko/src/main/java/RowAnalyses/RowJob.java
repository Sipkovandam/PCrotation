package RowAnalyses;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import Tools.FileUtils;

public abstract class RowJob
{
	private boolean hasSingleColHeader = false;
	String writeFn = null;
	String[] valueNames = null;
	
	public void execute(RowJobExecutor row, int lineNumber, int threadNumber)
	{
	}
		
	public String getWriteFn()
	{
		return writeFn;
	}

	public void setWriteFn(String writeFn)
	{
		this.writeFn = writeFn;
	}
	
	public boolean checkIncludeElement()
	{
		 return true;
	}

	public void writeLine(	int lineNumber,
											String writeLine, RowJobExecutor rowExecutor, boolean hasSingleColHeader, int threadNumber) throws IOException, InterruptedException
	{
		rowExecutor.write(writeLine, this.getWriteFn(), hasSingleColHeader, lineNumber, threadNumber);
	}

	public String[] getValueNames()
	{
		return valueNames;
	}

	public void setValueNames(String[] valueNames)
	{
		this.valueNames = valueNames;
	}

	protected void writeResult(double value, RowJobExecutor rowExecutor, int lineNumber, int threadNumber) throws FileNotFoundException, IOException, InterruptedException
	{
		String writeLine = rowExecutor.getRowName(threadNumber).concat("\t").concat(Double.toString(value).concat("\n"));
		writeLine(lineNumber, writeLine, rowExecutor, true, threadNumber);
	}
	
	protected void writeResult(int value, RowJobExecutor rowExecutor, int lineNumber, int threadNumber) throws FileNotFoundException, IOException, InterruptedException
	{
		String writeLine = rowExecutor.getRowName(threadNumber).concat("\t").concat(Integer.toString(value).concat("\n"));
		writeLine(lineNumber, writeLine, rowExecutor, true, threadNumber);
	}

	protected void writeResult(double[] values, RowJobExecutor rowExecutor, int lineNumber, int threadNumber) throws FileNotFoundException, IOException, InterruptedException
	{
		String writeLine = null;
		if(values!=null)
			writeLine = rowExecutor.getRowName(threadNumber).concat(FileUtils.doubleArrayToWriteString(values)).concat("\n");
		else
			writeLine="dontWriteMe";
		writeLine(lineNumber, writeLine, rowExecutor, false, threadNumber);		
	}

	public boolean hasSingleColHeader()
	{
		return hasSingleColHeader;
	}

	public void setHasSingleColHeader(boolean hasSingleColHeader)
	{
		this.hasSingleColHeader = hasSingleColHeader;
	}
}
