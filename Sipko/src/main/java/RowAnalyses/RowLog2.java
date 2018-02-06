package RowAnalyses;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;

import Tools.FileUtils;

public class RowLog2 extends RowJob
{
//	RowZscoreCalculator()
//	{
//		this.setValueNames(new String[]{zScoresName});
//	}
	
	public RowLog2()
	{
		this.setHasSingleColHeader(false);
	}
	
	@Override
	public void execute(RowJobExecutor rowExecutor, int lineNumber, int threadNumber)
	{
		try
		{		
			//log2 values
			double[] log2Values = log2(rowExecutor.getInputValues(threadNumber));
			super.writeResult(log2Values, rowExecutor, lineNumber, threadNumber);
		}catch(Exception e){e.printStackTrace();}
	}
	private double[] log2(double[] values) throws FileNotFoundException, IOException
	{
		double[] log2Values = new double[values.length];
		double logVal = Math.log(2);
		for(int v = 0; v < values.length; v++)
		{
			log2Values[v]= Math.log(values[v] + 1) / logVal;;
		}
		return log2Values;
	}

}
