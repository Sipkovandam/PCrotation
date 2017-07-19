package MatrixScripts;

import java.util.ArrayList;
import java.util.HashMap;
import java.io.BufferedReader;
import java.io.BufferedWriter;

import Tools.FileUtils;
import Tools.SampleSheetGafParser;
import Tools.Script;

public class LaneMerger extends Script<LaneMerger>
{
	//merges samples that were ran in separate lanes but that belong to the same run sample
	private String countsFnComment = "/root/folder/counts.txt; //File containing the counts for each sample in a different column. Some samples are run on multiple lanes and are merged in this script";
	private String countsFn = null;
	private String sampleSheetFnComment= "/root/folder/samplesheet.txt; //GAF sample sheet defining which samples where ran in which lanes";
	private String sampleSheetFn = null;
	private String laneIndexInFileNameComment = "4; //index the element containing the lane numbers after splitting on underscores";
	private Integer laneIndexInFileName = 4;//index the element containing the lane numbers after splitting based on "_"
	private String summedFnComment = "/root/folder/counts_LanesMerged.txt.gz; OPTIONAL //output file";
	private String summedFn = null;
	
	public void run()
	{
		try
		{
			if(summedFn==null)
				summedFn = FileUtils.removeExtention(countsFn) + "_LanesMerged.txt.gz";
			
			//parse GAF sample sheet
			SampleSheetGafParser gafSheet = new SampleSheetGafParser();
			gafSheet.setSampleSheetFn(sampleSheetFn);
			gafSheet.setSeparator(",");
			gafSheet.run();
			
			int nextIndex = 0;
			HashMap<Integer,Integer> oldToNewIndex = new HashMap<Integer,Integer>();
			HashMap<String,Integer> internalIdToNewIndex = new HashMap<String,Integer>();
			
			ArrayList<String> newNames = new ArrayList<String>();
			
			BufferedReader countsReader = FileUtils.createReader(countsFn); 
			
			String[] sampleNames = countsReader.readLine().split("\t");
			//fill table that contains conversion of the old column indexes to the new indexes
			for(int h = 1; h< sampleNames.length; h++)
			{
				int newIndex = 0;
				String internalId = gafSheet.sampleSheet.get("internalSampleID").get(sampleNames[h]);//get the internalId of this sample
				
				if(internalId!= null && internalIdToNewIndex.containsKey(internalId))
				{//if a sample is run on multiple lanes then
					newIndex=internalIdToNewIndex.get(internalId);
					String oldName = newNames.get(newIndex);
					
//					String[] oldInfo = oldName.split("_");
//					//String laneInfo = sampleNames[h].split("_")[laneIndexInFileName];
//					String newName = oldInfo[0];
//					
//					for(int l = 1; l < laneIndexInFileName; l++)
//						newName+="_"+oldInfo[l];
//					newName+="_";
//					//add lane element in the right place, you can do this if you want to keep all the lane numbers of the merged lanes in the filename
//					//newName = orderLaneNumbersAndAdd(newName,oldInfo,laneInfo);
					String newName =oldName.replaceAll("_L._", "_");
					
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
			
			BufferedWriter summedFileWriter = FileUtils.createWriter(summedFn);
			String line = null;
			for(String newName : newNames)
			{
				summedFileWriter.write("\t"+newName);
			}
			summedFileWriter.write("\n");
			
			while((line = countsReader.readLine())!=null)
			{
				double[] newValues = new double[nextIndex];
				String[] row = line.split("\t");
				for(int v = 1; v < row.length; v++)
				{
					newValues[oldToNewIndex.get(v)] += Double.parseDouble(row[v]);
				}
				
				String newLine = row[0];
				for(int v = 0; v < newValues.length; v++)
				{
					newLine=newLine.concat("\t");
					newLine=newLine.concat(String.valueOf(newValues[v]));
				}
				summedFileWriter.write(newLine.concat("\n"));
			}
			summedFileWriter.close();
			p("Summed file written to:" + summedFn);
		}catch(Exception e){e.printStackTrace();}
	}

	private String orderLaneNumbersAndAdd(String newName, String[] oldInfo, String laneInfo)
	{
		String[] lanes = oldInfo[laneIndexInFileName].split("(?=L)");
		boolean alreadyAdded = false;
		for(String lane : lanes)
		{
			if(lane.compareTo(laneInfo)<0 || alreadyAdded==true)
				newName += lane;
			else
			{
				newName += laneInfo+lane;
				alreadyAdded= true;
			}
		}
		if(!alreadyAdded)
			newName += laneInfo;
		
		newName +="_"+oldInfo[laneIndexInFileName+1];
		return newName;
	}
	
	public Integer getLaneIndexInFileName()
	{
		return laneIndexInFileName;
	}

	public void setLaneIndexInFileName(Integer laneIndexInFileName)
	{
		this.laneIndexInFileName = laneIndexInFileName;
	}
	
	public String getCountsFn()
	{
		return countsFn;
	}

	public void setCountsFn(String countsFn)
	{
		this.countsFn = countsFn;
	}

	public String getSampleSheetFn()
	{
		return sampleSheetFn;
	}

	public void setSampleSheetFn(String sampleSheetFn)
	{
		this.sampleSheetFn = sampleSheetFn;
	}

	public String getSummedFn()
	{
		return summedFn;
	}

	public void setSummedFn(String summedFn)
	{
		this.summedFn = summedFn;
	}
}
