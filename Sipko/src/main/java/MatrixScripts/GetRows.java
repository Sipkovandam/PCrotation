package MatrixScripts;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;

import Tools.FileUtils;
import Tools.Script;

public class GetRows extends Script<GetRows>
{
	String fileName = null;
	String rowsToGetFn = null;
	//String fileName2 = "ENSG00000268903,ENSG00000269981,ENSG00000225630";
	String writeName = null;
	String remove = null;
	boolean header = true;
	boolean originalOrder = false;
	
	@Override
	public void run()
	{
		try
		{
			init();
			
			MyMatrix rowsToInclude = getIncludeRowMatrix();
			
			BufferedReader reader = FileUtils.createReader(this.fileName);
			BufferedWriter writer = FileUtils.createWriter(this.writeName);
			String line = null;
			if(header)
			{
				line=reader.readLine();
				writer.write(line+"\n");//write the first line by default (headers)
			}
			Hashtable<String, Integer> toGet = rowsToInclude.namesToHash(rowsToInclude.rowNames);
			String[] results = new String[toGet.size()];
			
			int l = 0;
			while((line = reader.readLine()) != null)
			{	
				if(remove != null)
					line = line.replace(remove, "");
				String rowName = line.split("\t")[0];
				if(toGet.containsKey(rowName))
				{
					if(originalOrder)
						results[l] = line;
					else
						results[toGet.get(rowName)] = line;
					//writer.write(line+"\n");
					
					l++;
				}
			}
			
			for(int r = 0; r < results.length; r++)
			{
				if(results[r] == null)
					continue;
				writer.write(results[r]+"\n");
			}
			
			writer.close();
			reader.close();
			System.out.println("File written to:" + this.writeName);
		}catch(Exception e){e.printStackTrace();}
	}

	private void init()
	{
		if(this.writeName==null)
			this.writeName=FileUtils.addBeforeExtention(this.fileName,"_selectedRows");
		
	}

	private MyMatrix getIncludeRowMatrix()
	{
		MyMatrix file2=null;
		if(!this.rowsToGetFn.contains(",") && this.rowsToGetFn.contains(".txt"))
			file2 = new MyMatrix(this.rowsToGetFn,true,false);
		else{
			String[] rowNames = this.rowsToGetFn.split(",");
			file2 = new MyMatrix(rowNames.length, 1);
			file2.rowNames = rowNames;
			file2.colNames[0] = "-";
		}
		return file2;
	}

	public String getFileName()
	{
		return fileName;
	}

	public void setFileName(String fileName)
	{
		this.fileName = fileName;
	}

	public String getRowsToGetFn()
	{
		return rowsToGetFn;
	}

	public void setRowsToGetFn(String rowsToGetFn)
	{
		this.rowsToGetFn = rowsToGetFn;
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

	public boolean isOriginalOrder()
	{
		return originalOrder;
	}

	public void setOriginalOrder(boolean originalOrder)
	{
		this.originalOrder = originalOrder;
	}
}
