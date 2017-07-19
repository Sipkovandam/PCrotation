package Listeners.FileNameListeners;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.HashMap;

import PizzlyClasses.Fusion;
import PizzlyClasses.PizzlyFusionStructure;
import Tools.FileUtils;
import Tools.Script;

public class PizzlyFusionFileOperator extends FileNameListener
{
	HashMap<String,int[]> fusionCounts;
	
	@Override
	public void run(String pizzlyFusionFn)
	{
		try
		{
			PizzlyFusionStructure pizzlyFusionStructure = new PizzlyFusionStructure();
			pizzlyFusionStructure=(PizzlyFusionStructure) pizzlyFusionStructure.readVars(pizzlyFusionFn);
			BufferedWriter sampleFusionWriter = FileUtils.createWriter(FileUtils.removeExtention(pizzlyFusionFn)+ "_table.txt"); 
			
			for(Fusion fusion : pizzlyFusionStructure.getFusions())
			{
				String fusionName = getFusionName(fusion);
				addCountsToFusionCountsHash(fusion, fusionName);
				writeFusionToFile(sampleFusionWriter, fusion, fusionName);
			}
			sampleFusionWriter.close();
		}catch(Exception e){e.printStackTrace();}
	}

	private void writeFusionToFile(	BufferedWriter sampleFusionWriter,
									Fusion fusion, String fusionName) throws IOException
	{
		String overlappingReads = Integer.toString(fusion.getPaircount());
		sampleFusionWriter.write(fusionName.concat("\t").concat(overlappingReads).concat("\n"));
	}

	private String getFusionName(Fusion fusion)
	{
		String geneA = fusion.getGeneA().get("id");
		String geneB = fusion.getGeneB().get("id");
		
		String fusionName = geneA.concat("_").concat(geneB);
		return fusionName;
	}

	private void addCountsToFusionCountsHash(Fusion fusion, String fusionName)
	{		
		int[] counts = this.fusionCounts.get(fusionName);
		if(counts == null)
		{
			counts = new int[]{0,0};
		}
		//number of samples
		counts[0]++;
		//number of reads overlapping the splice junction
		counts[1]+= fusion.getPaircount();
		
		this.fusionCounts.put(fusionName, counts);
	}
	
	public HashMap<String, int[]> getFusionCounts()
	{
		return fusionCounts;
	}

	public void setFusionCounts(HashMap<String, int[]> fusionCounts)
	{
		this.fusionCounts = fusionCounts;
	}
}
