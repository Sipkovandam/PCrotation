package STAR;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashSet;

import MatrixScripts.GetRows;
import Tools.FileUtils;
import Tools.Script;

public class OutlierSpliceSiteRetriever extends Script<OutlierSpliceSiteRetriever>
{
	String outlierSpliceFolder = null;
	String rawMatrixFn = null;
	String correctedMatrixFn = null;
	int spliceNameCol = 1;
	
	@Override
	public void run()
	{
		try
			{
			File[] files = new File(outlierSpliceFolder).listFiles(); 
			HashSet<String> spliceSites = getRareSpliceSites(files); 
		
			String rowsToGet = convertToString(spliceSites);
			
			String resultFolder = FileUtils.makeFolderNameEndWithSlash(outlierSpliceFolder)+"AllOutlierSites/";
			new File(resultFolder).mkdir();
			
			//initate getrows from raw file
			GetRows getRowsRaw = new GetRows();
			getRowsRaw.setFileName(rawMatrixFn);
			getRowsRaw.setRowsToGetFn(rowsToGet.toString());
			getRowsRaw.setWriteName(resultFolder+new File(rawMatrixFn).getName());

			//initate getrows from corrected file
			GetRows getRowsCorrected = new GetRows();
			getRowsCorrected.setFileName(correctedMatrixFn);
			getRowsCorrected.setRowsToGetFn(rowsToGet.toString());
			getRowsCorrected.setWriteName(resultFolder+new File(correctedMatrixFn).getName());
			
			getRowsRaw.run();
			getRowsCorrected.run();
			
		}catch(Exception e){e.printStackTrace();}
	}

	private String convertToString(HashSet<String> spliceSites)
	{
		StringBuilder rowsToGet = new StringBuilder();
		boolean first = true;
		for(String spliceSite : spliceSites)
		{
			if(!first)
				rowsToGet.append(",");
			rowsToGet.append(spliceSite);
			first = false;
		}
		return rowsToGet.toString();
	}

	private HashSet<String> getRareSpliceSites(File[] files) throws FileNotFoundException, IOException
	{
		HashSet<String> spliceSites = new HashSet<String>();
		for(File file : files)
		{
			if(file.isDirectory() || !file.getName().contains(".txt"))
				continue;
			BufferedReader reader = FileUtils.createReader(file.getAbsolutePath());
			
			String line = reader.readLine(); //get rid of header
			while((line = reader.readLine())!= null)
			{
				String[] eles = line.split("\t");
				String spliceName = eles[spliceNameCol];
				spliceSites.add(spliceName);
			}
		}
		return spliceSites;
	}

	public String getInputFolder()
	{
		return outlierSpliceFolder;
	}

	public void setInputFolder(String inputFolder)
	{
		this.outlierSpliceFolder = inputFolder;
	}

	public String getOutlierSpliceFolder()
	{
		return outlierSpliceFolder;
	}

	public void setOutlierSpliceFolder(String outlierSpliceFolder)
	{
		this.outlierSpliceFolder = outlierSpliceFolder;
	}

	public String getRawMatrixFn()
	{
		return rawMatrixFn;
	}

	public void setRawMatrixFn(String rawMatrixFn)
	{
		this.rawMatrixFn = rawMatrixFn;
	}

	public String getCorrectedMatrixFn()
	{
		return correctedMatrixFn;
	}

	public void setCorrectedMatrixFn(String correctedMatrixFn)
	{
		this.correctedMatrixFn = correctedMatrixFn;
	}

	public int getSpliceNameCol()
	{
		return spliceNameCol;
	}

	public void setSpliceNameCol(int spliceNameCol)
	{
		this.spliceNameCol = spliceNameCol;
	}
}

