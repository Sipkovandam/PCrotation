package STAR;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;

import MatrixScripts.MyMatrix;
import Tools.FileUtils;
import Tools.Script;

public class RatioCalculator extends Script<RatioCalculator>
{
	//Calculates ratios based on expression (line by line)
	String expressionFolder = null;
	String writeFn = null;
	boolean spliceSitesOnRows = true;
	
	public void run()
	{
		
		try
		{
			if(writeFn == null)
				writeFn = new File(expressionFolder).getParent()+"/ratios.txt.gz";
			
			File folder = new File(expressionFolder);
			File[] expressionFiles = folder.listFiles();
			
			p("Files listed, calculating ratios and merging files");
			String header = "\t";
			
			int i = 0;
			for(File expressionFN : expressionFiles)
			{
				if(i%1000==0)
					System.out.println(i+"/"+ expressionFiles.length+ " done");
				if(expressionFN.isDirectory())
					continue;
				
				MyMatrix expression = new MyMatrix(expressionFN.getAbsolutePath());
				
				//create matrix with the ratios instead of expression values per splice for each sample
				//Splice sites are written on the rows, samples on the columns
				MyMatrix ratios = null;
				if(spliceSitesOnRows)
					ratios = expression.toRatiosPerCol(true);
				else
				{
					ratios = expression.toRatiosPerRow(true);
					ratios.transpose();
				}

				if(i==0)//write the file as is
					ratios.write(writeFn);
				else //append it onto the file
					ratios.write(writeFn, true);
				i++;
			}
		}catch(Exception e){e.printStackTrace();}
	}

	public String getExpressionFolder()
	{
		return expressionFolder;
	}

	public void setExpressionFolder(String expressionFolder)
	{
		this.expressionFolder = expressionFolder;
	}

	public String getWriteFn()
	{
		return writeFn;
	}

	public void setWriteFn(String writeFn)
	{
		this.writeFn = writeFn;
	}
}
