package Slurm;

public class SlurmVariables extends ClusterVariables
{
	private final String fileHeader = "#!/bin/bash";
	private final String jobName = "#SBATCH --job-name=";
	private final String logsFolder = "#SBATCH --output=";
	private final String errorsFolder = "#SBATCH --error=";
	private final String walltime = "#SBATCH --time=";
	private final String threads = "#SBATCH --threads=";
	private final String maxMemory = "#SBATCH --mem=";
	
	private final String extra ="#SBATCH --nodes=1\n"+
						"#SBATCH --export=NONE\n" + 
						"#SBATCH --get-user-env=L\n";
	
	private final String loadModule = "ml ";
	private final String loadModule2 = ""; 
	
	@Override
	public String getLoadModule()
	{
		return loadModule;
	}

	@Override
	public String getExtra()
	{
		return extra;
	}

	@Override
	public String getFileHeader()
	{
		return fileHeader;
	}

	@Override
	public String getJobName()
	{
		return jobName;
	}

	@Override
	public String getLogsFolder()
	{
		return logsFolder;
	}

	@Override
	public String getErrorsFolder()
	{
		return errorsFolder;
	}

	@Override
	public String getWalltime()
	{
		return walltime;
	}

	@Override
	public String getThreads()
	{
		return threads;
	}

	@Override
	public String getMaxMemory()
	{
		return maxMemory;
	}

	@Override
	public String getLoadModule2()
	{
		return loadModule2;
	}
}
