package PCA_testcases;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

import org.apache.commons.io.FileUtils;
import org.junit.Assert;

public class CompareFiles 
{
	public static void compare(String one, String two)
	{
		compare(one, two, true);
	}
	
	public static void compare(String one, String two, boolean ignoreFirstField)//only use this on small files! (I want it to print the whole thing if it is incorrect)
	{
		String first = null;
		String second = null;
		try {
			first = read(one);
			second = read(two);
		} catch (IOException e) {e.printStackTrace();}
		 
		if(ignoreFirstField && first.split("\t").length>1)
		{
			first=first.substring(first.indexOf("\t"), first.length());
			second=second.substring(second.indexOf("\t"), second.length());
		}
		
		//remove anything over 5 decimals.. (there are rounding error that I do not care about...)
		String regex = "(.[0-9]{5})[0-9]*([E\t\n])";
		String replace =  "$1$2";
		first = first.replaceAll(regex, replace);
		second = second.replaceAll(regex, replace);
		
//		first = first.trim();
//		second = second.trim();

//		System.out.println(first);	
//		System.out.println("\n");
//		System.out.println(second);
		Assert.assertEquals(first,second);
	}
	public static String read(String path) throws IOException
	{
		BufferedReader reader = null;
		
		if(path.endsWith(".gz"))
			reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(path))));
		else
			reader = new BufferedReader(new FileReader(new File(path)));

		String file = "";
		String line = null;
		while((line=reader.readLine())!=null)
			file+=line+"\n";
		reader.close();
		return file;
	}
	
}
