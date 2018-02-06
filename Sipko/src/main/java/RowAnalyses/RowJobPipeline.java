package RowAnalyses;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import Tools.FileUtils;
import Tools.Script;

public class RowJobPipeline extends Script<RowJobPipeline>
{
	private String comments = "//all variables in this file containing (Comment) are comments on the variable below and do not need to be initiated";
	String fileNameComment = "/root/directory/valueMatrix.txt; OPTIONAL; Input file for which the rowAverages should be calculated";
	String fileName = null;
	String writeFolder= null;
	String stepListComment = "//e.g. averages,center,removeNoVarianceRows,log2; MANDATORY //comma separated list of all steps to take in this order";
	String stepList = "";
	
	int nThreads=1;
	
	@Override
	public void run()
	{
		try
		{
			writeFolder=FileUtils.makeFolderNameEndWithSlash(writeFolder);
			FileUtils.makeDir(writeFolder);
			String[] steps = stepList.toLowerCase().replace(" ", "").split(",");
			String newName = new File(fileName).getName();
			
			String currentName = fileName;

			RowJobExecutor rowJobExecutor = new RowJobExecutor(this.nThreads);
			rowJobExecutor.setWriteFolder(FileUtils.makeFolderNameEndWithSlash(this.writeFolder));
			for(String step: steps)
			{
				
				String fileForNextStep = currentName;
				log("step = " + step);
				switch (step)
				{
					case "averages":
						RowAverageCalculator rowAverageCalculator = new RowAverageCalculator();
						rowAverageCalculator.setWriteFn(FileUtils.removeExtention(newName)+"_rowAverages.txt");
						rowJobExecutor.addJob(rowAverageCalculator);
						continue;
					case "center":
						RowCenterer center = new RowCenterer();
						newName=getNewName("_centered.txt.gz", currentName);
						center.setWriteFn(newName);
						rowJobExecutor.addJob(center);
						fileForNextStep=writeFolder+newName;
						break;
					case "removenovariancerows":
						RowNoVarianceRemover rowNoVarianceRemover = new RowNoVarianceRemover();
						newName=getNewName("_noVarRemoved.txt.gz", currentName);
						rowNoVarianceRemover.setWriteFn(newName);
						rowJobExecutor.addJob(rowNoVarianceRemover);
						fileForNextStep=writeFolder+newName;
						break;
					case "log2":
						RowLog2 rowLog2 = new RowLog2();
						newName=getNewName("_log2.txt.gz", currentName);
						rowLog2.setWriteFn(newName);
						rowJobExecutor.addJob(rowLog2);
						fileForNextStep=writeFolder+newName;
						break;
				}
				log("Executing " + rowJobExecutor.getRowJobs().size() + " jobs");
				RowJobExecutor.useExecutorsOnFile(rowJobExecutor, currentName);
				
				if(useNewFileNow(currentName,fileForNextStep))
				{
					currentName=fileForNextStep;
					rowJobExecutor = new RowJobExecutor(this.nThreads);
					rowJobExecutor.setWriteFolder(FileUtils.makeFolderNameEndWithSlash(this.writeFolder));
				}
			}
			
			
			
			System.out.println("Done! Files written to: " + writeFolder);
		}catch(Exception e){e.printStackTrace();}
		
	}

	private boolean useNewFileNow(String string)
	{
		// TODO Auto-generated method stub
		return false;
	}

	private boolean useNewFileNow(String currentName, String fileForNextStep)
	{
		if(currentName.equals(fileForNextStep))
			return false;
		return true;
	}

	private String getNewName(	String addString, String fileName)
	{
		String fn = new File(fileName).getName();
		log("fn=\t" + fn);
		String newName = FileUtils.removeExtention(fn)+ addString;
		return newName;
	}
}
