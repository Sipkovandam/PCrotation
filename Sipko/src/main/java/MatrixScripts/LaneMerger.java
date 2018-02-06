package MatrixScripts;

import java.util.ArrayList;
import java.util.HashMap;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;

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
			
			ArrayList<String> newNames = new ArrayList<String>();
			HashMap<Integer,Integer> oldToNewIndex = gafSheet.getOldToNewIndexes(countsFn, newNames);
			

			BufferedWriter summedFileWriter = FileUtils.createWriter(summedFn);
			BufferedReader countsReader = FileUtils.createReader(countsFn); 
			String line = null;
			for(String newName : newNames)
			{
				summedFileWriter.write("\t"+newName);
			}
			summedFileWriter.write("\n");
			
			while((line = countsReader.readLine())!=null)
			{
				double[] newValues = new double[newNames.size()];
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
			log("Summed file written to:" + summedFn);
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
