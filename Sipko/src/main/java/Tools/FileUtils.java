package Tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

public class FileUtils
{
	public static Hashtable<String,Double> readDoublehash(String fileName) throws FileNotFoundException, IOException
	{
		return readDoublehash(fileName, 0, 1);
	}
	public static Hashtable<String,Double> readDoublehash(String fileName, int col1, int col2) throws FileNotFoundException, IOException
	{
		Hashtable<String,Double> hash = new Hashtable<String,Double>();
		BufferedReader reader = createReader(fileName);
		String line = null;
		while((line = reader.readLine())!=null)
		{
			String[] cells = line.split("\t");
			hash.put(cells[col1], Double.parseDouble(cells[col2].replace(",", "")));
		}
		return hash;
	}
	
	public static ArrayList<String> readArrayList(String fileName) throws FileNotFoundException, IOException
	{
		return readArrayList(fileName, false);
	}
	public static ArrayList<String> readArrayList(String fileName, boolean firstColumnOnly) throws FileNotFoundException, IOException
	{
		
		ArrayList<String> lines = new ArrayList<String>();
	
		BufferedReader reader = null;
		
		if(fileName.endsWith(".gz"))		
			reader = new BufferedReader(new InputStreamReader(new GZIPInputStream (new FileInputStream(fileName))));
		else
			reader = new BufferedReader(new FileReader(new File(fileName)));
		
		String line = null;
		while((line=reader.readLine())!=null)
		{
			if(firstColumnOnly)
				line = line.split("\t")[0];
			lines.add(line);
		}
		return lines;
	}
	public static BufferedReader createReader(String fn) throws FileNotFoundException, IOException {
		BufferedReader reader = null;
		if(fn.endsWith(".gz"))
			reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fn))));
		else
			reader = new BufferedReader(new FileReader(new File(fn)));
		return reader;
	}
	
	public static BufferedWriter createWriter(String shellFN) throws FileNotFoundException, IOException {
		BufferedWriter writer = null;
		if(shellFN.endsWith(".gz"))
			writer = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(shellFN))));
		else
			writer = new BufferedWriter(new FileWriter(new File(shellFN)));
		return writer;
	}
	public static void makeDir(String folderName)
	{
		File folder = new File(folderName);
		if(!folder.exists())
			folder.mkdirs();
	}
	public static String replaceEnd(String string, String addition) {
		String result = string.replace(".txt", "").replace(".gz", "")+addition;
		return result;
	}
	public static void copy(String kallistoShellFN, String tempKallistoFN) throws FileNotFoundException, IOException {
		ClassLoader.getSystemClassLoader();
		
		ClassLoader loader = FileUtils.class.getClassLoader();
        System.out.println(loader.getResource("foo/FileUtils.class"));

		InputStream inputStream = ClassLoader.
		        getSystemResourceAsStream(kallistoShellFN);
		System.out.println("?=" + new File(kallistoShellFN).exists());
		System.out.println(inputStream);
		System.out.println(kallistoShellFN);
		InputStreamReader streamReader = new InputStreamReader(inputStream, "UTF-8");
		BufferedReader reader = new BufferedReader(streamReader);
		BufferedWriter writer = createWriter(tempKallistoFN);
		String line = null;
		while((line=reader.readLine())!=null)
		{
			System.out.println(line);
			writer.write(line+"\n");
		}
	}
	public static Hashtable<String, Integer> makeHash(String line) 
	{
		String[] cells = line.split("\t");
		Hashtable<String, Integer> index = new Hashtable<String, Integer>();
		for(int c = 0; c < cells.length; c++)
		{
			index.put(cells[c], c);
		}
		return index;
	}
}
