package MatrixScripts;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

import Tools.FileUtils;
import Tools.Script;

public class RowStatisticsGetter extends Script<RowStatisticsGetter>
{
	// calculates the averages and standard deviations per row, excluding NaN and INFINITE numbers.	
	String fileNameComment = "/root/directory/averges.txt; OPTIONAL; Input file for which the rowAverages should be calculated";
	String fileName = null;
	String writeFnComment = "/root/directory/averges.txt; OPTIONAL; Path where the file containing the row averages and standard deviations will be written";
	String writeFn = null;
	String absoluteComment = "false; OPTIONAL; Set to true if you wish to calculate the average of the absolute values instead";
	boolean absolute = false;
	transient int[] includeIndexes = null;

	public void run()
	{
		try
		{
			if(writeFn == null)
				writeFn = fileName.replace(".txt", "").replace(".gz", "")+"_rowAverages.txt";
			
			BufferedReader reader = FileUtils.createReader(this.fileName);
			BufferedWriter writerAverages = new BufferedWriter(new FileWriter(new File(this.writeFn)));
			
			String line = reader.readLine();//skip header row
			writerAverages.write("RowName\tAverage\tStandard deviation\n");
			
			if(this.includeIndexes ==null)
				this.includeIndexes=getAllColumnIndexesExceptFirst(line);
			
			while((line = reader.readLine()) !=null)
			{
				String[] eles = line.split("\t");
				double total = 0;
				double totalEles = (eles.length-1);//number of values
	
				int nValid = 0;//number of values that should be included in the average/stdev calculation
				double[] rowValues = new double[eles.length];//rowValues parsed from doubles; Indexes that should be excluded are set to NaN.
				for(int i = 0; i < this.includeIndexes.length; i++)// this makes this script a lot slower having to parse the same number 2X, but w/e
				{
					int e=this.includeIndexes[i];
					double number = Double.parseDouble(eles[e]);
					if(checkInclude(number,e))
					{
						rowValues[e]=number;
						nValid++;
					}
					else
					{
						rowValues[e]=Double.NaN;
					}
				}
				
				double[] rowValuesIncluded = new double[nValid];//used to calculate stdevs
				int index = 0;
				
				for(int i = 0; i < this.includeIndexes.length; i++)// this makes this script a lot slower having to parse the same number 2X, but w/e
				{
					int e=this.includeIndexes[i];
					double number = rowValues[e];
					if(!checkInclude(number,e))
					{
						totalEles--;
						continue;
					}
					if(this.absolute)
						total+=Math.abs(number);
					else
						total+=number;
					
					rowValuesIncluded[index] = number;
					index++;
				}
				double average = total/totalEles;
				double variance = org.apache.commons.math3.stat.StatUtils.variance(rowValuesIncluded);
				double standardDev = java.lang.Math.pow(variance,0.5);
				writerAverages.write(eles[0]+"\t" + average+"\t" + standardDev+ "\n");
			}
			writerAverages.close();
			reader.close();

			System.out.println("Done! File written to: " + writeFn);
		}catch(Exception e){e.printStackTrace();}
	}
	
	private int[] getAllColumnIndexesExceptFirst(String line)
	{
		String[] eles = line.split("\t");
		int[] includeIndexes=new int[eles.length-1];
		for(int e = 1;e< eles.length;e++)
			includeIndexes[e-1]=e;
		return includeIndexes;
	}

	private boolean checkInclude(double number, int e)
	{
		if(!Double.isNaN(number) && !Double.isInfinite(number))
			return true;
		return false;
	}

	public String getFileName()
	{
		return fileName;
	}
	public void setFileName(String fileName)
	{
		this.fileName = fileName;
	}
	public String getWriteFn()
	{
		return writeFn;
	}
	public void setWriteFn(String writeFn)
	{
		this.writeFn = writeFn;
	}
	public boolean isAbsolute()
	{
		return absolute;
	}
	public void setAbsolute(boolean absolute)
	{
		this.absolute = absolute;
	}

	public int[] getIncludeIndexes()
	{
		return includeIndexes;
	}

	public void setIncludeIndexes(int[] includeIndexes)
	{
		this.includeIndexes = includeIndexes;
	}
}
