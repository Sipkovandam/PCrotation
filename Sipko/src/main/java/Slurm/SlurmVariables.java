package Slurm;

public class SlurmVariables extends ClusterVariables
{
	String fileHeader = "#!/bin/bash";
	String jobName = "#SBATCH --job-name=";
	String logsFolder = "#SBATCH --output=";
	String errorsFolder = "#SBATCH --errorsFolder=";
	String walltime = "#SBATCH --walltime=";
	String threads = "#SBATCH --threads=";
	String maxMemory = "#SBATCH --maxMemory=";
	String nodes = "#SBATCH --nodes=";
	
	String extra =		"#SBATCH --export=NONE\n" + 
						"#SBATCH --get-user-env=L\n";
}
