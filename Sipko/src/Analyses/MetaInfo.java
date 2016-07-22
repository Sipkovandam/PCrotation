package Analyses;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.zip.GZIPInputStream;

public class MetaInfo
{
	//decided to do this in R instead; faster and can make nice pictures ;)
	public Hashtable<String, String[]> columns = new Hashtable<String, String[]>();
	public Hashtable<String, Integer> rowIndex = new Hashtable<String, Integer>();
	public Hashtable<Integer, String> colIndex = new Hashtable<Integer, String>();
	
	MetaInfo(String fn)
	{
		readFile(fn);//assumes the first cell also has a value (it skips this field)
	}
	private void readFile(String fn) 
	{
		try
		{
			int nLines = getNumberOfLines(fn);
			
			BufferedReader reader = getReader(fn);

			String line = null;
			int r = -1;
			while((line=reader.readLine())!=null)
			{
				String[] eles = line.split("\t");
				if(r==-1)//if it is the header column
				{
					for(int c = 1; c < eles.length; c++)
						colIndex.put(c,eles[c]);
					r++;
					continue;
				}
				
				rowIndex.put(eles[0], r);
				for(int c = 1; c < eles.length; c++)
				{
					String[] column = columns.get(c);
					if(column == null)
						column = new String[nLines];
					column[r]= eles[c];
					columns.put(colIndex.get(c),column);
				}
				r++;
			}
			
			reader.close();
		}catch(Exception e){e.printStackTrace();}
	}
	private int getNumberOfLines(String fn) throws IOException {
		BufferedReader reader = getReader(fn);
		int n = 0;
		while(reader.readLine() != null)
		{
			n++;
		}
		return n;
	}
	public BufferedReader getReader(String fn) throws FileNotFoundException, IOException
	{
		BufferedReader reader = null;
		if(fn.endsWith(".gz"))
		{
			GZIPInputStream inputStream = new GZIPInputStream(new FileInputStream(fn));
			reader = new BufferedReader(new InputStreamReader(inputStream,"US-ASCII"));
		}
		else
			reader = new BufferedReader(new FileReader(new File(fn)));
		return reader;
	}
	
}
