package Slurm;

public abstract class ClusterVariables
{
	private String fileHeader;
	private String jobName;
	private String logsFolder;
	private String errorsFolder;
	private String walltime;
	private String threads;
	private String maxMemory;
	private String loadModule;
	
	private String extra;
	private String moduleLoad;
	public String getExtra()
	{
		// TODO Auto-generated method stub
		return this.extra;
	}
	public String getFileHeader()
	{
		return fileHeader;
	}
	public void setFileHeader(String fileHeader)
	{
		this.fileHeader = fileHeader;
	}
	public String getJobName()
	{
		return jobName;
	}
	public void setJobName(String jobName)
	{
		this.jobName = jobName;
	}
	public String getLogsFolder()
	{
		return logsFolder;
	}
	public void setLogsFolder(String logsFolder)
	{
		this.logsFolder = logsFolder;
	}
	public String getErrorsFolder()
	{
		return errorsFolder;
	}
	public void setErrorsFolder(String errorsFolder)
	{
		this.errorsFolder = errorsFolder;
	}
	public String getWalltime()
	{
		return walltime;
	}
	public void setWalltime(String walltime)
	{
		this.walltime = walltime;
	}
	public String getThreads()
	{
		return threads;
	}
	public void setThreads(String threads)
	{
		this.threads = threads;
	}
	public String getMaxMemory()
	{
		return maxMemory;
	}
	public void setMaxMemory(String maxMemory)
	{
		this.maxMemory = maxMemory;
	}
	public String getLoadModule()
	{
		return loadModule;
	}
	public void setLoadModule(String loadModule)
	{
		this.loadModule = loadModule;
	}
	public String getModuleLoad()
	{
		return moduleLoad;
	}
	public void setModuleLoad(String moduleLoad)
	{
		this.moduleLoad = moduleLoad;
	}
	public void setExtra(String extra)
	{
		this.extra = extra;
	} 
}
