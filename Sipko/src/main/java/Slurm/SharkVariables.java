package Slurm;

public class SharkVariables extends ClusterVariables
{
	private final String fileHeader = "#!/bin/bash\n"+
						"#$ -S /bin/bash";
	private final String jobName = "#$ -N ";
	private final String logsFolder = "#$ -o ";
	private final String errorsFolder = "#$ -e ";
	private final String walltime = "#$ -l h_rt=";
	private final String threads = "#$ -pe BWA ";
	private final String maxMemory = "#$ -l h_vmem=";

	private final String extra = "#$ -q all.q\n";
	
	private final String loadModule = "export PATH=";
	private final String loadModule2 = ":$PATH"; 

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
