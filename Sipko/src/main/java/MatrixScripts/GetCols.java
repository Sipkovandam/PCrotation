package MatrixScripts;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;

import Tools.FileUtils;
import Tools.Script;

public class GetCols extends Script<GetCols>{
	//get some columns from a file works with really large files
	String fileName = null;
	String fileName2 = null;//column names that should be kept (should be a matrix with the names you want to keep on the rows). Alternatively a comma separated list of names can be supplied.
	String writeName = null;
	String remove = null;
	int extraColumns = 0;//keep the first few columns extra (columns after rowname)
	
	public void run()
	{
		try
		{
			init();
			
			//get names to keep
			MyMatrix namesToKeepOnRows = getNamesToKeep(this.fileName2);
			
			//Read the first line and get the new indexes of the columns that should be kept
			BufferedReader inputReader = FileUtils.createReader(fileName);
			String header = inputReader.readLine();
			String[] colNames= header.split("\t");

			BufferedWriter writer = FileUtils.createWriter(writeName);
			ArrayList<Integer> keepColNumbers = writeColumnHeaders(colNames,writer, namesToKeepOnRows, extraColumns);
			
			String line = null;		
			while((line=inputReader.readLine())!=null)
			{
				String[] rowCells= line.split("\t");
				//write the rowName
				writer.write(rowCells[0]);
				
				for(int col : keepColNumbers)
				{
					writer.write("\t");
					writer.write(rowCells[col]);
				}
				writer.write("\n");
			}
			log("Closing reader and writer");
			inputReader.close();
			writer.close();
			log("Done! File written to:" + writeName);
		} catch (FileNotFoundException e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private MyMatrix getNamesToKeep(String fileName22)
	{
		MyMatrix namesToKeepOnRows = null;
		if(!fileName2.contains(",") && fileName2.contains(".txt"))
			namesToKeepOnRows = new MyMatrix(this.fileName2,true,false);
		else{
			String[] rowNames = fileName2.split(",");
			namesToKeepOnRows = new MyMatrix(rowNames.length, 1);
			namesToKeepOnRows.rowNames = rowNames;
			namesToKeepOnRows.colNames[0] = "-";
		}
		return namesToKeepOnRows;
	}

	private void init()
	{
		if(writeName==null)
			writeName=FileUtils.removeExtention(fileName)+"_colsFrom_"+ FileUtils.removeExtention(new File(fileName2).getName()) +".txt.gz";
	}

	private ArrayList<Integer> writeColumnHeaders(String[] colNames, BufferedWriter writer, MyMatrix namesToKeepOnRows, int extraColumns2) throws IOException
	{
		ArrayList<Integer> keepColNumbers = new ArrayList<Integer>();
	
		for(int n = 1; n < colNames.length; n++)
		{
			if(remove!=null)
				colNames[n] = colNames[n].replaceAll(remove, "");
			
			log(colNames[n] +"\t"+ namesToKeepOnRows.getRowHash().containsKey((colNames[n])) );
			//add the extra few columns if there are any
			if(!namesToKeepOnRows.getRowHash().containsKey((colNames[n])) && n > extraColumns2)
				continue;
			writer.write("\t");
			writer.write(colNames[n]);
			keepColNumbers.add(n);
		}
		writer.write("\n");	
		return keepColNumbers;
	}

	public String getFileName()
	{
		return fileName;
	}

	public void setFileName(String fileName)
	{
		this.fileName = fileName;
	}

	public String getFileName2()
	{
		return fileName2;
	}

	public void setFileName2(String fileName2)
	{
		this.fileName2 = fileName2;
	}

	public String getWriteName()
	{
		return writeName;
	}

	public void setWriteName(String writeName)
	{
		this.writeName = writeName;
	}

	public String getRemove()
	{
		return remove;
	}

	public void setRemove(String remove)
	{
		this.remove = remove;
	}

	public int getExtraColumns()
	{
		return extraColumns;
	}

	public void setExtraColumns(int extraColumns)
	{
		this.extraColumns = extraColumns;
	}
}
