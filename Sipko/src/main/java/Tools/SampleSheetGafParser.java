package Tools;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class SampleSheetGafParser extends Script<SampleSheetGafParser>
{
	/**
	 * 
	 */
	private static final long serialVersionUID = -6648836498592663166L;
	public transient HashMap<String,HashMap<String,String>> sampleSheet = new HashMap<String,HashMap<String,String>>();//<gafColumnName,HashMap<sampleName,entry>
	private String sampleSheetFn = null;
	
	public void run()
	{
		try
		{
			BufferedReader sampleSheetReader = FileUtils.createReader(sampleSheetFn);
			HashMap<Integer,String> colNumberToColname= new HashMap<Integer, String>();
			HashMap<String,Integer> colNameToColNumber= new HashMap<String, Integer>();
			
			String header=sampleSheetReader.readLine();
			p("header = " +header);
			String[] colNames = header.split(separator);
			for(int c = 0; c < colNames.length; c++)
			{
				String colName = colNames[c];
				sampleSheet.put(colName, new HashMap<String,String>());
				colNumberToColname.put(c, colName);
				colNameToColNumber.put(colName, c);
			}
			
			String line = null;
			while((line=sampleSheetReader.readLine())!=null)
			{
				String[] dataPerCol=line.split(separator);
				int sequencingStartDateCol = colNameToColNumber.get("sequencingStartDate");
				int sequencerCol = colNameToColNumber.get("sequencer");
				int runCol = colNameToColNumber.get("run");
				int flowcellCol = colNameToColNumber.get("flowcell");
				int lane = colNameToColNumber.get("lane");
				int barcode = colNameToColNumber.get("barcode");
				
				String sampleName = dataPerCol[sequencingStartDateCol]+"_"+
						dataPerCol[sequencerCol]+"_"+
						addZeros(dataPerCol[runCol])+"_"+
						dataPerCol[flowcellCol]+"_L"+
						dataPerCol[lane]+"_"+
						dataPerCol[barcode];
				
				for(int c = 0; c < colNames.length; c++)
				{
					String colName = colNumberToColname.get(c);
					HashMap<String,String> colEntries = sampleSheet.get(colName);
					colEntries.put(sampleName, dataPerCol[c]);					
				}
			}
			
		}catch(Exception e)
		{
			e.printStackTrace();
		}
	}
	public void loadSampleSheet()
	{
		this.run();
	}

	private String addZeros(String run)
	{
		while(run.length()<3)
			run = "0"+run;
		return run;
	}
	public String getSampleSheetFn()
	{
		return sampleSheetFn;
	}
	public void setSampleSheetFn(String sampleSheetFn)
	{
		this.sampleSheetFn = sampleSheetFn;
	}
	public String getSeparator()
	{
		return separator;
	}
	public void setSeparator(String separator)
	{
		this.separator = separator;
	}

	private String separator=",";

	public HashMap<Integer, Integer> getOldToNewIndexes(String  countsFn, ArrayList<String> newNames) throws IOException
	{
		HashMap<String,Integer> internalIdToNewIndex = new HashMap<String,Integer>();
		HashMap<Integer,Integer> oldToNewIndex = new HashMap<Integer,Integer>();
		
		BufferedReader countsReader = FileUtils.createReader(countsFn); 
		
		String[] sampleNames = countsReader.readLine().split("\t");

		
		int nextIndex = 0;
		//fill table that contains conversion of the old column indexes to the new indexes
		for(int h = 1; h< sampleNames.length; h++)
		{
			int newIndex = 0;
			String internalId = this.sampleSheet.get("internalSampleID").get(sampleNames[h]);//get the internalId of this sample
			
			if(internalId!= null && internalIdToNewIndex.containsKey(internalId))
			{//if a sample is run on multiple lanes then
				newIndex=internalIdToNewIndex.get(internalId);
				String oldName = newNames.get(newIndex);
				
//				String[] oldInfo = oldName.split("_");
//				//String laneInfo = sampleNames[h].split("_")[laneIndexInFileName];
//				String newName = oldInfo[0];
//				
//				for(int l = 1; l < laneIndexInFileName; l++)
//					newName+="_"+oldInfo[l];
//				newName+="_";
//				//add lane element in the right place, you can do this if you want to keep all the lane numbers of the merged lanes in the filename
//				//newName = orderLaneNumbersAndAdd(newName,oldInfo,laneInfo);
				String newName =oldName.replaceAll("_L[0-9]{1,2}_", "_");
				
				newNames.set(newIndex, newName);
			}
			else
			{
				newNames.add(sampleNames[h]);
				newIndex = nextIndex;
				internalIdToNewIndex.put(internalId, newIndex);
				nextIndex++;
			}
			oldToNewIndex.put(h, newIndex);
		}
		return oldToNewIndex;
	}
	public HashMap<String, String> getCol(String string)
	{
		return sampleSheet.get(string);
	}
}
