package Tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import PCA.Var;

public class JSONutil <T>
{
 	
	public T read(String jsonFN, T object)//reads a jsonObject file into an object
	{
		T var = object;
		String jsonString = "";
		try
		{
			File file = new File(jsonFN);
			if(!file.exists())
			{
				System.out.println("Json file does not exist:" + jsonFN);
				System.exit(0);
			}
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String line = null;
			while((line=reader.readLine())!=null)
			{
				jsonString+=line;
			}
			reader.close();
		}catch(Exception e){e.printStackTrace();}
		
		Gson gson = new GsonBuilder().setPrettyPrinting().create();
		var = gson.fromJson(jsonString, object.getClass());
		System.out.println(gson.toJson(var));
		
		return var;
	}
	public void write(String jsonFN, T variables)
	{
		Gson gson = new GsonBuilder().setPrettyPrinting().create();
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
