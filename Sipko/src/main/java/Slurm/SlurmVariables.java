package Slurm;

public class SlurmVariables extends ClusterVariables
{
	private final String fileHeader = "#!/bin/bash";
	private final String jobName = "#SBATCH --job-name=";
	private final String logsFolder = "#SBATCH --output=";
	private final String errorsFolder = "#SBATCH --errorsFolder=";
	private final String walltime = "#SBATCH --walltime=";
	private final String threads = "#SBATCH --threads=";
	private final String maxMemory = "#SBATCH --maxMemory=";
	
	private final String extra ="#SBATCH --nodes=1\n"+
						"#SBATCH --export=NONE\n" + 
						"#SBATCH --get-user-env=L\n";
	
	private final String loadModule = "ml ";
	
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
}
