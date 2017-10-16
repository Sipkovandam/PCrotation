package Tools;

import java.io.File;
import java.lang.reflect.Constructor;
import java.lang.reflect.Method;
import Tools.Script;

public class Run
{
	public static void main(String[] args) throws Exception
	{
		run(args);
	}

	public static void run(String[] args) throws Exception
	{
		if (args.length != 1 && args.length != 3)
		{
			System.out.println("This splice site pipeline requires a json input file. Options:\n" + "1. Supply the pathname to the json input file, for example:\n" + "</root/directory/config.json>\n" + "2. Request for an empty configuration file for a certian script:\n " + "getjson package.classname </root/directory/config.json> \n");
			System.exit(1);
		}

		if (args.length == 3)//create config file
		{
			createConfigFile(args);
		}
		if (args.length == 1)//create config file
		{
			executeScript(args[0]);
		}
	}

	public static void executeScript(String jsonFN) throws Exception
	{
		File file = new File(jsonFN);
		if (!file.exists())
		{
			System.out.println("Config file does not exist: " + jsonFN);
			System.out.println("Exiting");
			System.exit(1);
		}
		if (file.isDirectory())
		{
			System.out.println("Please supply a  filename, not a directory");
			System.out.println("Exiting");
			System.exit(1);
		}

		String classToRun = getClassToRun(jsonFN);
		System.out.println("Class to run:\t" + classToRun);
		Class<Script<?>> runClass = (Class<Script<?>>) Class.forName(classToRun);//add if runnable and endable
		Constructor<Script<?>> runClassConstructor = runClass.getDeclaredConstructor();
		Script<?> objectLoader = runClassConstructor.newInstance();
		Method readVars = runClass.getMethod(	"readVars",
												String.class,
												runClass.getClass());

		Script<?> runObject = (Script<?>) readVars.invoke(	objectLoader,
															jsonFN,
															runClass);
		Method startMethod = runClass.getMethod("start",												String.class);
		startMethod.invoke(	runObject,
							jsonFN);
		Method runMethod = runClass.getMethod("run");
		runMethod.invoke(runObject);
		Method endMethod = runClass.getMethod("end");
		endMethod.invoke(runObject);
	}

	public static String getClassToRun(String jsonFN)
	{
		Script<?> runClassJson = new Script<Object>().read(jsonFN);
		return runClassJson.getClassName();
	}

	public static void createConfigFile(String[] args) throws Exception
	{
		File dir = new File(args[2]);
		File parentDir = new File(dir.getParent());
		System.out.println(parentDir.getAbsolutePath());
		if (!parentDir.exists())
		{
			System.out.println("Parent directory does not exist: " + parentDir.getAbsolutePath() + "\n" + "Exiting");
			System.exit(1);
		}
		if (!dir.exists())
		{
			System.out.println("Directory does not exist: " + dir.getAbsolutePath() + "\n" + "Creating directory");
			dir.mkdir();
		}
		
		
		System.out.println("Creating config file");
		String classToRun = args[1];
		System.out.println("ExecutorClass=" + classToRun);
		//get the class
		Class<Script<?>> configClass = (Class<Script<?>>) Class.forName(classToRun);
		//get the constructor
		Constructor<Script<?>> constructor = configClass.getDeclaredConstructor();
		System.out.println(constructor);
		//create an object using the constructor;
		Script<?> configWriter = constructor.newInstance();
		//get the method
		Method writeConfig = configClass.getMethod(	"createConfig",
													String.class);
		//invoke the method
		writeConfig.invoke(	configWriter,
							dir.getAbsolutePath());
	}
}
