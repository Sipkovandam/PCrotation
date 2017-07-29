package RowAnalyses;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import Tools.FileUtils;

public abstract class RowJob
{
	String writeFn = null;
	String[] valueNames = null;
	
	public void execute(RowJobExecutor row, int lineNumber)
	{
	}
	
	public void executeOnInitiation(RowJobExecutor row, double value)
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

	public BufferedWriter getWriter(RowJobExecutor rowExecutor, boolean hasSingleColheader) throws FileNotFoundException, IOException
	{	
		HashMap<String, BufferedWriter> rowExecutorFileWriters = rowExecutor.getFileWriters();
		BufferedWriter fileWriter = rowExecutorFileWriters.get(this.getWriteFn());
		
		if(fileWriter == null)
		{
			String colheader = null;
			if(hasSingleColheader)
				colheader ="\t"+FileUtils.removeExtention(getWriteFn());
			else
			{
				colheader = FileUtils.StringArrayToWriteString(rowExecutor.getDataColHeaders());
			}
			
			String fileName = rowExecutor.getWriteFolder()+this.getWriteFn();
			fileWriter = FileUtils.createWriter(fileName);
			String header = "Rowname"+colheader+"\n";
			fileWriter.write(header);
			rowExecutorFileWriters.put(this.getWriteFn(), fileWriter);
		}
		return fileWriter;
	}

//	public void writeResultSlots(RowJobExecutor rowExecutor, boolean hasSingleColHeader) throws IOException
//	{
//		String line = null;
//		while(resultSlots.size() > this.resultsWrittenIndex && (line = resultSlots.get(this.resultsWrittenIndex))!=null)//get the next line that should be written. If it is available yet, write it and move to the next line to write
//		{
//			resultSlots.set(this.resultsWrittenIndex, null);
//			this.resultsWrittenIndex++;
//			
//			BufferedWriter writer = getWriter(rowExecutor, hasSingleColHeader);
//			writer.write(line);	
//		}
//	}
//
	public void writeLine(	int lineNumber,
											String writeLine, RowJobExecutor rowExecutor, boolean hasSingleColHeader) throws IOException
	{
		BufferedWriter writer = getWriter(rowExecutor,hasSingleColHeader);
		writer.write(writeLine);		
	}

	public String[] getValueNames()
	{
		return valueNames;
	}

	public void setValueNames(String[] valueNames)
	{
		this.valueNames = valueNames;
	}
}
