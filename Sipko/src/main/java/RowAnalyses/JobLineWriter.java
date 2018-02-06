package RowAnalyses;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.concurrent.atomic.AtomicBoolean;

import Tools.FileUtils;

public class JobLineWriter
{
	String writeFolder = null;
	String writeFn = null;
	
	int currentRow = 0;
	BufferedWriter writer = null;
	final int writeBufferSize = 1000;//should not be smaller than number of threads.// if this is too small it could cause the buffer to fill so quickly the variable is not set by the time it reiterates over the buffer (Just don't make it smaller than 1000, no piont)
	String[] linesToWrite = new String[writeBufferSize];

	boolean isWriting= false;
	private final Object lock = new Object();

	public JobLineWriter(	String writeFolder,
							String writeFn,
							String[] dataColHeaders, boolean hasSingleColHeader) throws FileNotFoundException, IOException
	{
		this.writeFn = writeFn;
		this.writeFolder = writeFolder;
		
		this.writer = getWriter(this.writeFolder, this.writeFn, dataColHeaders, hasSingleColHeader);
		
		// TODO Auto-generated constructor stub
	}

	public BufferedWriter getWriter(String writeFolder, String writeFn, String[] dataColHeaders, boolean hasSingleColHeader) throws FileNotFoundException, IOException
	{	
		BufferedWriter fileWriter = writer;
		
		if(fileWriter == null)
		{
			String colheader = null;
			if(hasSingleColHeader || dataColHeaders ==null)
				colheader ="\t"+FileUtils.removeExtention(writeFn);
			else
			{
				colheader = FileUtils.StringArrayToWriteString(dataColHeaders);
			}
			
			String fileName = FileUtils.makeFolderNameEndWithSlash(writeFolder)+writeFn;
			fileWriter = FileUtils.createWriter(fileName);
			String header = "Rowname"+colheader+"\n";
			fileWriter.write(header);
		}
		return fileWriter;
	}

	public void close() throws IOException
	{
		if(writer != null)
			writer.close();
	}

	public void writeIntoBufferAndWriteBuffer(	int lineNumber,
								String writeLine, int threadNumber) throws IOException, InterruptedException
	{
		int bufferPosition = lineNumber%writeBufferSize;
		
		while(linesToWrite[bufferPosition] !=null)//buffer is full
		{
			Thread.sleep(100);
		}
		
		linesToWrite[bufferPosition]=writeLine;
//		System.out.println("lineNumber= " + lineNumber + "\tccurrentrow = " + currentRow+ "\t"+  linesToWrite[lineNumber] );
	}

	public void writeLinesInBuffer() throws IOException, InterruptedException
	{
		String line = linesToWrite[currentRow];
		while(line != null)
		{
			linesToWrite[currentRow] = null;
			if(!line.equals("dontWriteMe"))//value in case a line is not to be written
			{
				writer.write(line);
			}
			currentRow++;
			if(currentRow>=writeBufferSize)
				currentRow = 0;
			line = linesToWrite[currentRow];
		}
	}
	
	
}


