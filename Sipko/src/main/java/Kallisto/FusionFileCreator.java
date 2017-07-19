package Kallisto;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import Listeners.FileNameListeners.FileNameListener;
import Listeners.FileNameListeners.PizzlyFusionFileMaker;
import Listeners.FileNameListeners.PizzlyFusionFileOperator;
import Tools.FileSearcher;
import Tools.FileUtils;
import Tools.Script;

public class FusionFileCreator extends Script<FusionFileCreator>
{
	FileSearcher kallistoFusionFileSearcher = new FileSearcher();
	
	PizzlyFusionFileMaker pizzlyExecutor = new PizzlyFusionFileMaker();
	
	String writeFn = null;
	
	@Override
	public void run()
	{
		try
		{
			//find all files
			kallistoFusionFileSearcher.run();
			
			//loop over files and run them through pizzly
			String fusionFilesFn = kallistoFusionFileSearcher.getWriteName();
			List<FileNameListener> jobList = new ArrayList<FileNameListener>();
			jobList.add(pizzlyExecutor);
			//executeJobsOnFiles(fusionFilesFn, jobList);
			
			p("Pizzly Files created, searching output files");
			//Search the pizzly output files
			FileSearcher pizzlyFileSearcher = new FileSearcher();
			pizzlyFileSearcher.setFolders(kallistoFusionFileSearcher.getFolders());
			pizzlyFileSearcher.setSearchStrings(new String[]{"pizzly.json"});
			pizzlyFileSearcher.setWriteName(new File(kallistoFusionFileSearcher.getWriteName()).getParent()+"/pizzlyFileNames.txt");
			pizzlyFileSearcher.run();
			
			p("Mering pizzly splice files into 1 summary file");
			//Count the fusions in all the pizzly output files
			PizzlyFusionFileOperator pizzlyFusionFileOperator = new PizzlyFusionFileOperator();
			HashMap<String,int[]> fusionCounts = new HashMap<String,int[]>();
			pizzlyFusionFileOperator.setFusionCounts(fusionCounts);
			jobList = new ArrayList<FileNameListener>();
			jobList.add(pizzlyFusionFileOperator);
			
			executeJobsOnFiles(pizzlyFileSearcher.getWriteName(), jobList);
			
			//write the resulting counts
			if(writeFn == null)
				writeFn = new File(fusionFilesFn).getParent()+"fusionSummary.txt";
			
			FileUtils.writeHashtable(fusionCounts, writeFn, "GeneFusion\tObservedInSamples\tTotalReadsSupportingFusion");
			p("Finished. File written to: \t"+ fusionCounts);
			
		}catch(Exception e){e.printStackTrace();}
	}

	private void executeJobsOnFiles(String fusionFilesFn,
									List<FileNameListener> jobList) throws IOException
	{
		BufferedReader fusionFilesFnReader = FileUtils.createReader(fusionFilesFn);
		fusionFilesFnReader.lines().forEach(fileName -> runJobs(fileName, jobList));
	}

	private void runJobs(	String fileName,
								List<FileNameListener> jobList)
	{
		for(FileNameListener fileNameJob : jobList)
		{
			fileNameJob.run(fileName);
		}
	}
}
