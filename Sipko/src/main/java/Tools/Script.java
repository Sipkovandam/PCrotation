package Tools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.lang.reflect.Modifier;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.Format;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.lang.SerializationUtils;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import PCA.RlogLargeMatrix;
import STAR.SpliceMerger_InfiniteFileSizes;
import STAR.STAR_Pipeline;
import Tests.TestObject;
import Tests.TestObject3;

public class Script <T> implements Tools.Runnable, Serializable
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	static transient final Format timeFormat = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
	
	
	private String className = this.getClass().getName();
	private String startTime = timeFormat.format(new Date());
	private String executor = Script.class.getProtectionDomain().getCodeSource().getLocation().toString();
	private transient double start = System.nanoTime();
	private String runTime = "Not finished";
	public String jsonFN = "No need to initiate this variable, but used in pipelines to identify the scripts that where used in all substeps";
	
	protected Script()
	{
		
	}
	
	public void writeConfig()
	{
		writeConfig(jsonFN, (T) this);
	}
	public void writeConfig(String jsonFN)
	{
		writeConfig(jsonFN, (T) this);
	}
	public void createConfig(String jsonFN, T object, boolean writeTransient, boolean serializeNulls) 
	{
		this.setJsonFN(jsonFN);
		if(new File(this.jsonFN).exists())
		{
			p("Json file already exists at: " + this.jsonFN +"\n"
			+ "Exiting");
			System.exit(1);
		}
		writeConfig(this.jsonFN, object,writeTransient,serializeNulls);	
	}
	public Script<T> read(String jsonFilename)//function should be in separate interface...
	{
		return readVars(jsonFilename, null);
	}
	public Script<T> readVars(String jsonFilename, Class<?> classType)//function should be in separate interface...
	{
		String jsonString = "";
		try
		{
			File file = new File(jsonFilename);
			if(!file.exists())
			{
				p("Json file does not exist:" + jsonFilename);
				System.exit(0);
			}
			jsonString = new String(Files.readAllBytes(Paths.get(jsonFilename)));
		}catch(Exception e){e.printStackTrace();}
		
		GsonBuilder gsonBuilder = new GsonBuilder().setPrettyPrinting(); 
		
		Gson gson = gsonBuilder.create();
		if(classType==null)
			classType=this.getClass();
		Script<T> var = gson.fromJson(jsonString, classType);
		var.jsonFN=jsonFilename;
		return var;
	}
		
	public List<Object> r(Object ... objects)//retuns a list of objects... Handy for passing objects from a function
	{
		List<Object> objectList = new ArrayList<>();
		for(Object object : objects)
			objectList.add(object);
		return objectList;
	}
	
	public Script<T> clone(String jsonFN) throws CloneNotSupportedException
	{
		Script<T> clone = (Script<T>) SerializationUtils.clone(this);//deep clone
		this.jsonFN=jsonFN;
		clone.jsonFN=jsonFN;
		writeConfig();
		return clone;
	}
	
	public void run()
	{
		p("run() not defined for this class:" + className);
		p("exiting");
		System.exit(2);
	}
	public void run(String[] args) throws Exception
	{
		p("run(String[] args) not defined for this class:" + className);
	}
	
	public void writeConfig(String jsonFN, T script)
	{
		writeConfig(jsonFN, script,true,true);	
	}
	public void writeConfig(String jsonFN, boolean writeTransient, boolean serializeNulls)
	{
		writeConfig(jsonFN,(T) this, writeTransient, serializeNulls);
	}
	public void createConfig(String jsonFN)
	{
		createConfig(jsonFN, (T) this,true,true);	
	}
	public String getNewJsonFN() {
		return getNewJsonFN(null);
	}
	
	public String getNewJsonFN(String addString) {
		if(addString == null)
			addString="";
		else
			addString+="_";
		String jsonFN = "config_"+getClassName()+ "_" +addString+getDate()+".json";
		return jsonFN;
	}
	public String getJsonFN() {
		return jsonFN;
	}
	public void setJsonFN(String jsonFN) 
	{
		File file = new File(jsonFN);
		if(file.isDirectory())
			if(!file.exists())
				p("Directory does not exist, exiting: " + file);
			else
				this.jsonFN=FileUtils.makeFolderNameEndWithSlash(jsonFN)+this.getNewJsonFN(null);
		else
			this.jsonFN=jsonFN;
	}

	public void writeConfig(String jsonFN, T script, boolean writeTransient, boolean serializeNulls)
	{
//		if(jsonFN.startsWith("null"))
//			return;
		if(this.jsonFN.contains("No need to initiate this variable, but used in pipelines to identify the scripts that where used in all substeps"))
			this.jsonFN=getNewJsonFN();
		this.setJsonFN(jsonFN);
		
		this.executor = Script.class.getProtectionDomain().getCodeSource().getLocation().toString();
		
		GsonBuilder gsonBuilder = new GsonBuilder().setPrettyPrinting();
		if(serializeNulls)
			gsonBuilder = gsonBuilder.serializeNulls();
		
		Gson gson = gsonBuilder.create();
		//p(gson.toJson(this));
		//p(gson.toString());
		
		try
		{
			BufferedWriter writer = new BufferedWriter(new FileWriter(new File(this.getJsonFN())));
			
			write_StripNestedScriptVariables(gson.toJson(script), writer);
			writer.close();
		}catch(Exception e){e.printStackTrace();}
		p("Config file written at:\t"+this.getJsonFN());
	}
	
	private void write_StripNestedScriptVariables(String jsonString, BufferedWriter writer) throws IOException {
		HashMap<String, Integer> stringAllowence = new HashMap<String,Integer>();
				
		//write inverse
		ArrayList<String> strippedJsonLines=this.stripLines(jsonString,getStringAllowence());
		
		//strippedJsonLines are in reverse order
		for(int l = strippedJsonLines.size()-1; l >= 0;l--)
		{
			writer.write(strippedJsonLines.get(l)+"\n");
		}
	}
	public HashMap<String, Integer> getStringAllowence() {//function to be overwritten by child classes
		HashMap<String, Integer> stringAllowence = new HashMap<String,Integer>();//number of times each variable can appear in final jsonFile
		stringAllowence.put("className", 1);
		stringAllowence.put("startTime", 1);
		stringAllowence.put("executor", 1);
		stringAllowence.put("runTime", 1);
		return stringAllowence;
	}

	private ArrayList<String> stripLines(String jsonString, HashMap<String, Integer> stringAllowence) {
		String[] lines = jsonString.split("\n");
		ArrayList<String> strippedJsonLines = new ArrayList<String>();
		//strip lines from bottom to top (so the outter class variables are still written)
		boolean stripComma =false;
		for(int l=lines.length-1; l>=0;l--)
		{
			String line = lines[l];
			if(stripComma && line.endsWith(","))
				line = line.substring(0, line.length()-1);
			
			if(!line.contains("\""))
			{
				strippedJsonLines.add(line);
				if(lines[l].endsWith("},") || lines[l].endsWith("}"))
					stripComma=true;
				continue;
			}
			
			String variable = line.split("\"")[1];
			int left = -1;
			if(stringAllowence.containsKey(variable))
				left = stringAllowence.get(variable);
			
			if(left!=0)
			{
				strippedJsonLines.add(line);
				left--;
				if(left>=0)
					stringAllowence.put(variable,left);
				stripComma=false;
			}		
		}
		return strippedJsonLines;
	}
	public void start()
	{
		start(null);
	}
	public void start(String jsonFN)
	{
		startTime = timeFormat.format(new Date());
		start = System.nanoTime();
		if(!jsonFN.contains("No need to initiate this variable, but used in pipelines to identify the scripts that where used in all substeps"))
		{
			this.jsonFN=jsonFN;
			writeConfig(this.jsonFN, true,true);
		}
		
//		if(this.jsonFN!=null)
//		{
//			this.setJsonFN(FileUtils.addBeforeExtention(this.jsonFN, "_usedConfigs"));
//			writeConfig(this.jsonFN, true,true);
//		}		
//		else
//			p("jsonFN not set, config file not written");
	}
	public void end()
	{
		this.runTime = this.getRunTime();
		writeConfig(this.jsonFN, true,true);
		p("Done! Runtime:\t" + this.runTime);
	}
	
	public String getRunTime()
	{
		double end = System.nanoTime();
		int miliSeconds = (int)((end-start)/1000/1000);
		int seconds = miliSeconds/1000;
		int minutes = seconds/60;
		int hours = minutes/60;
		int days = hours/24;
		
		String runTime = miliSeconds%1000 + " ms";
		if(seconds>0)
			runTime = seconds%60+" seconds "+runTime;
		if(minutes>0)
			runTime = minutes%60+" minutes "+runTime;
		if(hours>0)
			runTime = hours%60+" hours "+runTime;
		if(days>0)
			runTime = days%24+" days "+runTime;
		
		return runTime;
	}
	
	public void p(Object line)
	{
		String time = timeFormat.format(new Date());
        System.out.println(time + " (" + this.getClassName() + "):\t" + line);
	}

	public String getClassName() {
		if(this.className == null)
			this.className=this.getClass().getName();
		return className;
	}

	public static String getDate()
	{
		Format dateFormat = new SimpleDateFormat("yyyy-MM-dd");
		String date = dateFormat.format(new Date());
		return date;
	}
}

