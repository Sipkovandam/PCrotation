package MatrixScripts;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;

import Tools.FileUtils;
import Tools.Script;

public class CountValuesPerRow extends Script<CountValuesPerRow>
{
	String fnComment = "/root/directory/file.txt; Tab separated file from which the values should be counted per row";
	String fn = null;
	String writeFnComment = "/root/directory/file.txt; File where the results should be written";
	String writeFn = "";
	String valueComment = "0; Value that should be counted";
	int testVal = 0;
	String matrixStartColComment = "0; Column from which the matrix stars (the column with the row names)";
	int matrixStartCol = 0;
	
	public void run()
	{
		try
		{
			BufferedReader reader = FileUtils.createReader(fn);
			String line = reader.readLine();//skip header
			BufferedWriter writer = FileUtils.createWriter(writeFn);

			//write header
			writeHeader(testVal, line, writer, matrixStartCol);
				
			countValues(testVal, writer, reader, matrixStartCol);
			
			reader.close();
			writer.close();
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
	}

	private void countValues(	int testVal,
								BufferedWriter writer, BufferedReader reader, int matrixStartCol) throws NumberFormatException, IOException
	{
		String line = null;
		while((line=reader.readLine())!=null)
		{		
			String rowElements[] = line.split("\t");
			//write the rowname
			writer.write(rowElements[matrixStartCol]);
			String writeLine="";
			int tmpTestValue=testVal;
			
			double[] values = rowToValues(rowElements, matrixStartCol);
			
			while(true)
			{
				int n = 0;
				for(int e = 0; e < values.length; e++)
				{
					if(values[e] < tmpTestValue)
						n++;
				}
				writeLine=writeLine.concat("\t"+n);
				if(tmpTestValue<1)
					break;
				tmpTestValue/=2;
			}
			writeLine=writeLine.concat("\n");
			writer.write(writeLine);
		}
	}

	private double[] rowToValues(String[] rowElements, int matrixStartCol)
	{
		double[] values= new double[rowElements.length-1-matrixStartCol];
		for(int e = matrixStartCol+1; e < rowElements.length; e++)//1 because we want to skip the rowname
		{
			Double value = Double.parseDouble(rowElements[e]);
			values[e-1-matrixStartCol]=value;
		}
		return values;
	}

	private void writeHeader(int testVal, String line, BufferedWriter writer, int startCol2) throws IOException
	{
		String writeLine = "Rownames";
		while(true)
		{
			//write the rowname
			
			writeLine=writeLine.concat("\tNumber of values below "+testVal+" (max="+ (line.split("\t").length-1)+ ")");
			if(testVal<1)
				break;
			testVal/=2;
		}
		writeLine=writeLine.concat("\n");
		writer.write(writeLine);
	}
}
