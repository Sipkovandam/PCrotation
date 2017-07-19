package Listeners.FileNameListeners;

import java.io.File;

import Listeners.Listener;
import Tools.ExecCommand;
import Tools.FileUtils;

public class PizzlyFusionFileMaker extends FileNameListener
{
	private String pizzlyExecute = null;
	private String gtfFn = null;
	private String fastaFn = null;
	private String outputPostFix = null;// e.g. "pizzly" (file of interest is then named pizzly.json)
	private String pizzlyExtraArguments = null;

	@Override
	public void run(String fileName)
	{
		FileUtils.checkNull(this.pizzlyExecute, "pizzlyExecute");
		String outputDir = new File(fileName).getParent()+"/PizzlyOutput/";
		new File(outputDir).mkdir();
		checkArguments();
		
		//build the command line command
		String command = pizzlyExecute;
		command=command.concat(" -k 31 --cache index.cache.txt");
		command=command.concat(" --gtf ").concat(this.gtfFn);
		command=command.concat(" --fasta ").concat(fastaFn);
		command=command.concat(" --output ").concat(outputDir+outputPostFix);
		if(pizzlyExtraArguments!=null)
			command=command.concat(" ").concat(pizzlyExtraArguments);
		command=command.concat(" ").concat(fileName);
		
		ExecCommand exec = new ExecCommand(command);
		// p("execute output: \n"+exec.getOutput());
		if(exec.getError().length()>1)
		{
			p("Warning: Pizzly execute error: \t" + exec.getError());
			p("Warning: Command leading to error: \t" + command);
		}
		else
			p("Pizzly completed running with no errors on file:" + fileName);
	}

	public String getPizzlyExecute()
	{
		return pizzlyExecute;
	}

	public void setPizzlyExecute(String pizzlyExecute)
	{
		this.pizzlyExecute = pizzlyExecute;
	}

	public String getPizzlyArguments()
	{
		return pizzlyExtraArguments;
	}

	public void setPizzlyExtraArguments(String pizzlyArguments)
	{
		this.pizzlyExtraArguments = pizzlyArguments;
	}

	private void checkArguments()
	{
		if(pizzlyExecute ==null)
			System.out.println("Error <gtfFn> not set");
		if(gtfFn ==null)
			System.out.println("Error <gtfFn> not set");
		if(outputPostFix ==null)
			System.out.println("Error <outputPrefix> not set");
		
		if(gtfFn ==null || outputPostFix ==null || pizzlyExecute==null)
		{
			System.out.println("Exiting");
			System.exit(2);
		}
	}

	public String getGtfFn()
	{
		return gtfFn;
	}

	public void setGtfFn(String gtfFn)
	{
		this.gtfFn = gtfFn;
	}
	
	public String getPizzlyExtraArguments()
	{
		return pizzlyExtraArguments;
	}

	public String getOutputPrefix()
	{
		return outputPostFix;
	}

	public void setOutputPrefix(String outputPrefix)
	{
		this.outputPostFix = outputPrefix;
	}
}
