package RowAnalyses;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.apache.commons.math3.stat.descriptive.rank.Median;

import Tools.FileUtils;

public class RowNoVarianceRemover extends RowJob
{
	
	public RowNoVarianceRemover()
	{
		this.setHasSingleColHeader(false);
	}
	
	@Override
	public void execute(RowJobExecutor rowExecutor, int lineNumber, int threadNumber)
	{
		try
		{
			double[] values = rowExecutor.getInputValues(threadNumber);
			boolean hasVariance=false;
			for(double value : values)
			{
				if(value!=values[0])
				{
					hasVariance=true;
					break;
				}
			}

			if(hasVariance)
				super.writeResult(values, rowExecutor, lineNumber, threadNumber);
			else
				super.writeResult(null, rowExecutor, lineNumber, threadNumber);
		}catch(Exception e){e.printStackTrace();}
	}
}
