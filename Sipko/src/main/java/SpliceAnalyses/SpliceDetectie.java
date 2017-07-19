package SpliceAnalyses;

import java.util.ArrayList;
import java.util.List;

import MatrixScripts.MergeFiles;
import Tools.FileUtils;
import Tools.Runnable;
import Tools.Script;

public class SpliceDetectie extends Script<SpliceDetectie>
{
	String fn = null;
	String ensgGeneToLengthFn = null;
	
	public void run()
	{
		try
		{
			List<Runnable> steps = init();
			for(Runnable step : steps)
			{
				step.run();
			}
		}catch(Exception e){e.printStackTrace();}
	}

	private List<Runnable> init()
	{
		//initiate steps
		List<Runnable> steps = new ArrayList<>();
		//merge geneLengths
		MergeFiles fileMerger = new MergeFiles();
		fileMerger.setFn1(fn);
		fileMerger.setFn2(ensgGeneToLengthFn);
		fileMerger.setWriteFn(FileUtils.addBeforeExtention(fn, "_merged"));
		
		//get the count for the most used splice junction for each splice junction
		
		
		//add all steps
		steps.add(fileMerger);
		return steps;
	}
	
}
