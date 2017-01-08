package Tools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.text.Format;
import java.text.SimpleDateFormat;
import java.util.Date;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import PCA.RlogLargeMatrix;

public abstract class Script <T>
{
	protected String jsonFN = null;
	private final String className = this.getClass().getCanonicalName();//handy in the json file to see which script generated it
	private final String date = new SimpleDateFormat("dd/MM/yyyy HH:mm:ss").format(new Date());
	
	static transient final Format timeFormat = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
	private transient final double start = System.nanoTime();
	private String runTime = "Not finished";
	
	private void writeVars()
	{
		write(jsonFN, this);
	}
	
	public void write(String jsonFN, Script<T> script)
	{
		Gson gson = new GsonBuilder().serializeNulls().setPrettyPrinting().create();
		//System.out.println(gson.toJson(this));
		//System.out.println(gson.toString());
		System.out.println(jsonFN);
		try
		{
			BufferedWriter writer = new BufferedWriter(new FileWriter(new File(jsonFN)));
			
			writer.write(gson.toJson(script));
			writer.close();
		}catch(Exception e){e.printStackTrace();}
	}
	
	protected void end()
	{
		this.runTime = this.getRunTime();
		System.out.println("Runtime:\t" + this.runTime);
		writeVars();
	}
	
	private String getRunTime()
	{
		double end = System.nanoTime();
		int miliSeconds = (int)(end-start)/1000/1000;
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
	
	static protected void p(String line)
	{
		String time = timeFormat.format(new Date());
        System.out.println(time + "\t" + line);
	}
}
