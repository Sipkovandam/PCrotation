package Listeners.FileNameListeners;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.HashMap;

import PizzlyClasses.Fusion;
import PizzlyClasses.PizzlyFusionStructure;
import Tools.FileUtils;
import Tools.Script;

public class PizzlyFusionSummaryCreator extends FileNameListener
{
	HashMap<String,HashMap<String,Integer>> fusionCountsPerSample;//<fusion,<fusionCount,count>

	
	@Override
	public void run(String pizzlyFusionFn)
	{
		try
		{
			PizzlyFusionStructure pizzlyFusionStructure = new PizzlyFusionStructure();//contains all the fusions
			pizzlyFusionStructure=(PizzlyFusionStructure) pizzlyFusionStructure.read(pizzlyFusionFn);

			BufferedWriter fusionTableWriter = FileUtils.createWriter(FileUtils.removeExtention(pizzlyFusionFn)+ "_table.txt");
			
			int nIncludeIntable = 0;
			for(Fusion fusion : pizzlyFusionStructure.getFusions())
			{
				String fusionName = getFusionName(fusion);
				addCountsToFusionCountsHash(fusion, fusionName, pizzlyFusionFn);
			}
			log("Fusions:" + fusionCountsPerSample.size());

		}catch(Exception e){e.printStackTrace();}
	}



	private String getFusionName(Fusion fusion)
	{
		String geneA = fusion.getGeneA().get("id");
		String geneB = fusion.getGeneB().get("id");
		
		String fusionName = geneA.concat("_").concat(geneB);
		return fusionName;
	}

	private void addCountsToFusionCountsHash(Fusion fusion, String fusionName, String sampleName)
	{		
		HashMap<String,Integer> fusionCounts = this.fusionCountsPerSample.get(fusionName);
		if(fusionCounts==null)
			fusionCounts = new HashMap<String,Integer>();
		//number of reads overlapping the splice junction
		int counts = fusion.getPaircount();
		fusionCounts.put(sampleName, counts);
		
		this.fusionCountsPerSample.put(fusionName, fusionCounts);
	}
	
	public HashMap<String,HashMap<String,Integer>> getFusionCounts()
	{
		return fusionCountsPerSample;
	}

	public void setFusionCounts(HashMap<String,HashMap<String,Integer>> fusionCounts)
	{
		this.fusionCountsPerSample = fusionCounts;
	}
}
