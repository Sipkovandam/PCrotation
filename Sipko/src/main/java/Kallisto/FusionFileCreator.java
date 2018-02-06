package Kallisto;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import Listeners.FileNameListeners.FileNameListener;
import Listeners.FileNameListeners.PizzlyFusionFileMaker;
import Listeners.FileNameListeners.PizzlyFusionSummaryCreator;
import PizzlyClasses.Fusion;
import Tools.FileSearcher;
import Tools.FileUtils;
import Tools.Script;

public class FusionFileCreator extends Script<FusionFileCreator>
{
	FileSearcher kallistoFusionFileSearcher = new FileSearcher();
	
	PizzlyFusionFileMaker pizzlyExecutor = new PizzlyFusionFileMaker();
	
	String writeFn = null;
	int minReadSupportForInclusionInTable = 8;
	
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
			
			log("Pizzly Files created, searching output files");
			//Search the pizzly output files
			FileSearcher pizzlyFileSearcher = new FileSearcher();
			pizzlyFileSearcher.setFolders(kallistoFusionFileSearcher.getFolders());
			pizzlyFileSearcher.setSearchStrings(new String[]{"pizzly.json"});
			pizzlyFileSearcher.setWriteName(new File(kallistoFusionFileSearcher.getWriteName()).getParent()+"/pizzlyFileNames.txt");
			pizzlyFileSearcher.run();
			
			log("Mering pizzly splice files into 1 summary file");
			//Count the fusions in all the pizzly output files
			PizzlyFusionSummaryCreator pizzlyFusionFileOperator = new PizzlyFusionSummaryCreator();
			HashMap<String,HashMap<String,Integer>> fusionCounts = new HashMap<String,HashMap<String,Integer>>();//<fusion,<fusionCount,count>
			pizzlyFusionFileOperator.setFusionCounts(fusionCounts);
			
			jobList = new ArrayList<FileNameListener>();
			jobList.add(pizzlyFusionFileOperator);
			
			executeJobsOnFiles(pizzlyFileSearcher.getWriteName(), jobList);
			
			//write the resulting counts
			if(writeFn == null)
				writeFn = new File(fusionFilesFn).getParent()+"fusionSummary.txt";
			
			//FileUtils.writeHashtable(fusionCounts, writeFn, "GeneFusion\tObservedInSamples\tTotalReadsSupportingFusion");
			writeFusionToFile(fusionCounts, pizzlyFileSearcher.getWriteName());
			log("Finished. File written to: \t"+ pizzlyFileSearcher.getWriteName());
			
		}catch(Exception e){e.printStackTrace();}
	}
	
	private void writeFusionToFile(HashMap<String,HashMap<String,Integer>> fusionCounts, String sampleFileNames) throws IOException
	{
		BufferedWriter fusionSummaryWriter = FileUtils.createWriter(writeFn); 
		
		HashMap<String,Integer> fileNameToColumn = writeHeader(sampleFileNames, fusionSummaryWriter);
		
		log("Number of detected fusions:" + fusionCounts.keySet().size());
		for(String fusion : fusionCounts.keySet())
		{
			HashMap<String,Integer> sampleToCounts = fusionCounts.get(fusion);
			int[] counts = new int[fileNameToColumn.size()];
			int sum = 0;
			boolean include = false;
			for(String sample : sampleToCounts.keySet())
			{
				int count = sampleToCounts.get(sample);
				int index = fileNameToColumn.get(sample);
				counts[index]=count;
				sum+=count;
				
				if(count>0 && sample.matches(".*_0667_.*|.*_0694_.*|.*_0713_.*|.*_0096_.*"))
					include = true;
			}
			
			StringBuilder countBuilder = new StringBuilder();
			for(int count: counts)
			{
				countBuilder.append("\t"+count);
			}
			//if(sum >= minReadSupportForInclusionInTable)
			if(include)
				fusionSummaryWriter.write(fusion.concat(countBuilder.toString()).concat("\n"));	
		}
		
		
		fusionSummaryWriter.close();
	}

	private HashMap<String,Integer> writeHeader(String sampleFileNames, BufferedWriter fusionSummaryWriter)
	{
		HashMap<String,Integer> fileNameToColumn = new HashMap<String,Integer>();
		try
		{
			BufferedReader fusionFilesFnReader = FileUtils.createReader(sampleFileNames);
			String line = null;
			String header = "";
			
			int col = 0;
			while((line=fusionFilesFnReader.readLine())!=null)
			{
				String fileName = new File(new File(new File(line).getParent()).getParent()).getName();
				header+="\t"+fileName;
				fileNameToColumn.put(line, col);
				col++;
			}
			fusionSummaryWriter.write(header+"\n");
		}catch(Exception e){}
		return fileNameToColumn;
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
			if(new File(fileName).isDirectory())
				continue;
			fileNameJob.run(fileName);
		}
	}
}
