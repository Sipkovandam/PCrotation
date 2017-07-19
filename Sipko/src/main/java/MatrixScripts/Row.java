package MatrixScripts;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.descriptive.rank.Median;

import Tools.FileUtils;

public class Row
{
	public String rowName = null;
	public double[] values = null;
	Double avg = null;
	Double stdev = null;
	
	public double getAvg()
	{
		if(avg==null)
		{
			double[] stats = getStatistics();
			this.avg = stats[0];
			this.stdev = stats[1];
		}
		return avg;
	}

	public void setAvg(double avg)
	{
		this.avg = avg;
	}

	public double getStdev()
	{
		if(stdev==null)
		{
			double[] stats = getStatistics();
			this.avg = stats[0];
			this.stdev = stats[1];
		}
		return stdev;
	}

	public void setStdev(double stdev)
	{
		this.stdev = stdev;
	}

	public Row(String rowName, double[] values)
	{
		this.rowName=rowName;
		this.values = values;
	}
	
	public Row(	String rowName,
				int length)
	{
		this.rowName=rowName;
		this.values=new double[length];
	}
	
	public Row convertToZscores(Row row, double average, double stdev)
	{
		Row rowZscores = new Row(rowName, row.values.length);
		for(int v = 0; v < row.values.length; v++)
		{
			double zscore = (row.values[v]-average)/stdev;
			rowZscores.values[v] = zscore;
		}
		return rowZscores;
	}

	public static Row readRow(String line)
	{
		String[] eles = line.split("\t",2);
		
		String[] valuesString = eles[1].split("\t");
		double[] values = new double[valuesString.length];
		for(int v = 0; v < valuesString.length; v++)
		{
			values[v] = Double.parseDouble(valuesString[v]);
		}
		Row row = new Row(eles[0],values);
		row.rowName = eles[0];
		return row;
	}

	public String getRowName()
	{
		return rowName;
	}

	public void setRowName(String rowName)
	{
		this.rowName = rowName;
	}

	public double[] getValues()
	{
		return values;
	}

	public void setValues(double[] values)
	{
		this.values = values;
	}

	public double[] getStatistics()
	{
		return getStatistics(false);
	}
	
	public double[] getStatistics(boolean absolute)
	{
		return getStatistics(null,false);
	}
	
	public double[] getStatistics(int[] includeIndexes, boolean absolute)
	{
		int nValid = 0;//number of values that should be included int he calculation
		double total = 0;
		
		for(int i = 0; i < includeIndexes.length; i++)// this makes this script a lot slower having to parse the same number 2X, but w/e
		{
			int e = includeIndexes[i];
			double value = values[e];
			if(checkInclude(value,e))
			{
				nValid++;
			}
		}
		
		double[] rowValuesIncluded = new double[nValid];//used to calculate stdevs
		int index = 0;
		
		double totalEles = (rowValuesIncluded.length);//number of values
		//put the values to calculate the averages and stdevs over into 1 vector
		for(int i = 0; i < includeIndexes.length; i++)// this makes this script a lot slower having to parse the same number 2X, but w/e
		{
			int e = includeIndexes[i];
			double number = values[e];
			if(!checkInclude(number,e))
			{
				totalEles--;
				continue;
			}
			if(absolute)
				total+=Math.abs(number);
			else
				total+=number;

			rowValuesIncluded[index] = number;
			index++;
		}
		double average = total/totalEles;
		double variance = org.apache.commons.math3.stat.StatUtils.variance(rowValuesIncluded);
		double standardDev = java.lang.Math.pow(variance,0.5);
		double largest = org.apache.commons.math3.stat.StatUtils.max(rowValuesIncluded);
		double median = getMedian(rowValuesIncluded);
		return new double[]{average,standardDev,largest,median};
	}
	
	private boolean checkInclude(double number, int e)
	{
		if(!Double.isNaN(number) && !Double.isInfinite(number))
			return true;
		return false;
	}

	public void write(BufferedWriter zscoreWriter) throws IOException
	{
		zscoreWriter.write(this.rowName);
		for(int i = 0; i < values.length; i++)
		{
			zscoreWriter.write("\t");
			zscoreWriter.write(Double.toString(this.values[i]));
		}
		zscoreWriter.write("\n");
	}
	public double getMedian(double[] values)
	{
		Median median = new Median();
		double medianValue = median.evaluate(values);
		return medianValue;
	}

}
