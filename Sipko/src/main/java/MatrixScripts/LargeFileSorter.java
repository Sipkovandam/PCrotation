package MatrixScripts;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;

import Tools.FileUtils;
import Tools.Script;

public class LargeFileSorter extends Script<LargeFileSorter>
{
	private String fn;
	private String sortOrderFn;
	private String fnSorted;
	
	private String splitRegexRowNameAndTake2nd;

	@Override
	public void run()
	{
		try
		{
			if(fnSorted == null)
				fnSorted = FileUtils.addBeforeExtention(fn,"_sorted")+".txt";
			
			BufferedWriter sortedWriter = FileUtils.createWriter(fnSorted);
			writeHeader(fn, sortedWriter);
			
			LineSearcher fnLineSearcher = new LineSearcher(fn);
			
			BufferedReader sortedFnReader = FileUtils.createReader(sortOrderFn);
			
			String line = sortedFnReader.readLine();//skipheader
			long startTime = System.nanoTime();
			int l = 0;
			int lr = 0;
			while((line = sortedFnReader.readLine())!=null)
			{
				lr++;
				if(l%10000==0)
				{
					Long currentRunTime = (System.nanoTime()-startTime)/1000/1000/1000;
					System.out.println("lines Indexed = " + l + " \truntime=" + currentRunTime+ "seconds");
				}
				String rowName = line.split("\t",2)[0];
				
				String geneName = null;
				if(splitRegexRowNameAndTake2nd!=null)
				{
					String[] gene_Splice = rowName.split(splitRegexRowNameAndTake2nd);
					geneName = gene_Splice[0];
					rowName = gene_Splice[1];
					
					if(line.contains("__-"))
					{
						System.out.println("line failed=\t" + line);
						System.out.println("rowName=\t" + rowName);
						System.out.println("geneName=\t" + geneName);
						
						rowName=rowName.split("\\.")[1];
					}
				}

				String writeLine = fnLineSearcher.getLine(rowName);
				
				if(writeLine== null)
				{
					System.out.println("line failed=\t" + line);
					System.out.println("rowName=\t" + rowName);

					continue;
				}
				writeLine +="\n";
			
				if(splitRegexRowNameAndTake2nd!=null)
					writeLine = geneName+splitRegexRowNameAndTake2nd+writeLine;

				sortedWriter.write(writeLine);
				l++;
			}
			sortedWriter.close();
			System.out.println("Total lines read   = "+ lr);
			System.out.println("Total lines written   = "+ l);
		}catch(Exception e){e.printStackTrace();}		
	}


	private void writeHeader(	String fn,
								BufferedWriter sortedWriter) throws FileNotFoundException, IOException
	{
		String header = FileUtils.createReader(fn).readLine();
		sortedWriter .write(header+"\n");
	}


	public String getFn()
	{
		return fn;
	}


	public void setFn(String fn)
	{
		this.fn = fn;
	}


	public String getSortOrderFn()
	{
		return sortOrderFn;
	}


	public void setSortOrderFn(String sortOrderFn)
	{
		this.sortOrderFn = sortOrderFn;
	}


	public String getFnSorted()
	{
		return fnSorted;
	}


	public void setFnSorted(String fnSorted)
	{
		this.fnSorted = fnSorted;
	}
	
}
