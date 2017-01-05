package Tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import PCA.RlogLargeMatrix;
import PCA.Var;

public class JSONutil <T>
{
	public T read(String[] args, T var)
	{
		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase())
			{
				case "json":
					return read(value, var);
			}
		}
		return var;
	}
	
	public T read(String jsonFN, T var)//reads a jsonObject file into an object
	{
		String jsonString = "";
		try
		{
			File file = new File(jsonFN);
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String line = null;
			while((line=reader.readLine())!=null)
			{
				jsonString+=line;
			}
			reader.close();
		}catch(Exception e)
		{
			File file = new File(jsonFN);
			if(!file.exists())
				System.out.println("Json file does not exist:" + jsonFN);
			e.printStackTrace();
		}
		
		Gson gson = new GsonBuilder().setPrettyPrinting().create();
		var = gson.fromJson(jsonString, var.getClass());
		System.out.println(gson.toJson(var));
		
		return var;
	}
	public void write(String jsonFN, T variables)
	{
		Gson gson = new GsonBuilder().serializeNulls().setPrettyPrinting().create();
		//System.out.println(gson.toJson(this));
		//System.out.println(gson.toString());
		System.out.println(jsonFN);

		try
		{
			BufferedWriter writer = new BufferedWriter(new FileWriter(new File(jsonFN)));
			
			writer.write(gson.toJson(variables));
			writer.close();
		}catch(Exception e){e.printStackTrace();}
	}

	public String getFolderName(String fn) 
	{
		if(!fn.endsWith("/") && !fn.endsWith("\\"))
			fn = fn+"/";
		return fn;
	}
}
