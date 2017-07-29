package MatrixScripts;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.RandomAccessFile;
import java.nio.charset.StandardCharsets;
import java.util.HashMap;

import Tools.FileUtils;

public class LineSearcher
{
	String fileName = null;
	HashMap<String,Long> rowNameToFilePostion = null;
	RandomAccessFile randomAccesPointer = null;
	BufferedReader reader = null;
	BufferedInputStream bis =null;

	LineSearcher (String fileName)
	{
		this.fileName=fileName;
	}
	
	private HashMap<String,Long> loadlineNumberToFilePosition() throws FileNotFoundException, IOException
	{
		if(this.rowNameToFilePostion ==null)
		{
			String indexFn = fileName+".idx";
			if(!new File(indexFn).exists())
				createIndexFile(this.fileName, indexFn);
			this.rowNameToFilePostion=FileUtils.readStringToLongHash(indexFn);
		}
		return this.rowNameToFilePostion;
	}

	private void createIndexFile(	String fileName,
									String indexFn) throws IOException
	{
		if(fileName.endsWith(".gz"))
			System.out.println("File must be unzipped: " + fileName);
		randomAccesPointer = new RandomAccessFile(fileName, "r");
		FileInputStream fis = new FileInputStream(randomAccesPointer.getFD());
		BufferedInputStream bis = new BufferedInputStream(fis);
		reader = new BufferedReader( new InputStreamReader(bis, StandardCharsets.UTF_8));
		
		BufferedWriter indexWriter = FileUtils.createWriter(indexFn);
		
		long pos = randomAccesPointer.getFilePointer();
		
		StringBuilder writeLine = new StringBuilder();
		int l = 0;
		long startTime = System.nanoTime();
		String line = null;
		while((line=reader.readLine())!=null)
		{
			if(l%10000==0)
			{
				Long currentRunTime = (System.nanoTime()-startTime)/1000/1000/1000;
				System.out.println("lines Indexed = " + l + " \truntime=" + currentRunTime+ "seconds");
			}
			writeLine.setLength(0);
			String rowName = line.split("\t",2)[0];
			writeLine.append(rowName).append("\t").append(pos).append("\n");
			indexWriter.write(writeLine.toString());
			pos += line.length()+1;
			l++;		
		}
		indexWriter.close();
	}

	public String getLine(String rowName) throws FileNotFoundException, IOException
	{
		if(rowNameToFilePostion==null)
			rowNameToFilePostion=loadlineNumberToFilePosition();
		
		Long filePosition = rowNameToFilePostion.get(rowName);
		
		if(filePosition == null)
			return null;
		

		if(randomAccesPointer == null || bis ==null)
		{
			randomAccesPointer = new RandomAccessFile(fileName, "r");
			FileInputStream fis = new FileInputStream(randomAccesPointer.getFD());
			bis = new BufferedInputStream(fis);
			reader = new BufferedReader( new InputStreamReader(bis, StandardCharsets.UTF_8));
		}
		reader = new BufferedReader( new InputStreamReader(bis, StandardCharsets.UTF_8));
		randomAccesPointer.seek(filePosition);
		String line = reader.readLine();
		
		return line;
	}
}
