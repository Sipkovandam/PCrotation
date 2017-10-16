package Tests;

import java.io.BufferedReader;
import java.io.File;
import java.net.URI;
import java.nio.file.FileSystem;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.security.CodeSource;
import java.util.Collections;

import Tools.ExecCommand;
import Tools.FileUtils;

public class ReadIncludedFile
{
	
	public static void main(String[] args)
	{
		try
		{
			CodeSource src = ReadIncludedFile.class.getProtectionDomain().getCodeSource();
			System.out.println("src =" + src.getLocation().getPath());
			
			String kallistoShellFN = "resources/Kallisto200.sh";
			extractFile(src.getLocation().getPath(), kallistoShellFN);
			
			if(!new File(kallistoShellFN).exists())
				System.out.println("Slurmscript does not exist, please copy to:\t" + kallistoShellFN);
			
			if(new File(kallistoShellFN).exists())
				System.out.println("file exists: " + new File(kallistoShellFN).getAbsolutePath());
			
			
			

			
		}catch(Exception e){e.printStackTrace();}
	}

	private static void extractFile(String jarFn, String fileToGet)
	{
		String command = "jar xf " + jarFn + " " + fileToGet;
		System.out.println("Command =\t" + command);
		ExecCommand exec = new ExecCommand(command);
		// p("execute output: \n"+exec.getOutput());
		if(exec.getError().length()>1)
			System.out.println("execute error: \n" + exec.getError());
		else
			System.out.println("No errors occurred running the slurm scripts");
	}
}
