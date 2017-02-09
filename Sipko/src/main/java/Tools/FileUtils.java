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
import java.util.HashMap;
import java.util.Hashtable;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import java.net.URL;
import java.net.URLDecoder;

import org.apache.commons.io.IOUtils;

import PCA.RlogLargeMatrix;

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
		if(fileName == null || !new File(fileName).exists())
		{
			throw new FileNotFoundException("File does not exist: "+fileName);
		}
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
		String result = removeExtention(string)+addition;
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
	public static HashMap<String, String> readStringStringHash(String ensgToGeneSymbolFN) 
	{
		HashMap<String, String> conversionHash = new HashMap<String, String>();
		try 
		{
			BufferedReader reader = createReader(ensgToGeneSymbolFN);
			reader.lines().forEach(line -> 
			{
				String[] cells = line.split("\t");
				if(cells.length>1)
					conversionHash.put(cells[0], cells[1]);
			});
		} catch (FileNotFoundException e) {	e.printStackTrace();
		} catch (IOException e) {e.printStackTrace();}
		return conversionHash;
	}
	public static HashMap<String, Double> makeHash(String readCountFile, int i) {
		HashMap<String, Double> keyToValue = new HashMap<String, Double>();		
		try {
			BufferedReader reader = createReader(readCountFile);
			reader.lines().forEach(line -> 
			{
				String[] eles = line.split("\t");
				keyToValue.put(eles[0], Double.parseDouble(eles[i]));
			});
			reader.close();
		} catch (IOException e) {e.printStackTrace();}
		return keyToValue;
	}
	
	public static InputStream getResourceFile(String file) throws IOException
	  {
	    return Util.class.getResourceAsStream("/" + file);
	  }
	  
	  public static String getResourcePath(String file) throws IOException
	  {
	  	URL url =  Util.class.getResource("/" + file);
	  	if(url == null)
	  		throw new IOException("resource file not found: " + file);
	  	return url.toString().replaceAll("file:/","file:///");
	  }
	  
	  public static String readFile(String file) throws IOException
	  {
	    return IOUtils.toString(Util.getResourceFile(file));
	  }
	  
	  public static String prepareFileFromJar(String file) throws IOException
	  {
	    return FileUtils.prepareFileIntFromJar(file,false,true);
	  }
	  
	  public static String prepareBinaryFromJar(String file) throws IOException
	  {
	    return FileUtils.prepareFileIntFromJar(file,true,true);
	  }
	  
	  public static String getTempFile(String file) throws IOException
	  {
	    return FileUtils.prepareFileIntFromJar(file,false,false);
	  }
	  
	  public static void cleanTemp() throws IOException
	  {
	    String tempPath = FileUtils.getTempPath();
	    (new File(tempPath)).delete();
	  }
	  
	  private static String prepareFileIntFromJar(String file, boolean isBinary, boolean doCreate) throws IOException
	  {
	    // Copies to location next to JAR file
	    String tempPath = FileUtils.getTempPath();
	    (new File(tempPath)).mkdir();
	    File outBinFileName = new File(tempPath + file.substring(Math.max(0,file.lastIndexOf("/"))));
	    if (doCreate)
	    {
	      String fn = file.replaceAll("\\{OS\\}",getOs());
	      System.out.println(fn+"\t"+ new File(fn).exists() + "\t" + outBinFileName);
	      
	      InputStream stream = Util.getResourceFile(fn);
	      FileOutputStream output = new FileOutputStream(outBinFileName);
	      IOUtils.copy(stream, output);
	      output.close();
	      if (isBinary)
	        outBinFileName.setExecutable(true);
	      System.out.println("Creating file: " + outBinFileName.toString());
	    }
	    return outBinFileName.toString();
	  }
	  
	  private static String getTempPath() throws IOException
	  {
	    String path = Util.class.getProtectionDomain().getCodeSource().getLocation().getPath();
	    String jarPath = URLDecoder.decode(path,"UTF-8").toString();
	    return jarPath.substring(0,jarPath.lastIndexOf("/")) + "/temp/";
	  }

	  public static String getBigResource(Class clazz,String basePathName,String resourse) throws IOException
	  {
	    String path = clazz.getProtectionDomain().getCodeSource().getLocation().getPath();
	    String jarPath = URLDecoder.decode(path,"UTF-8").toString();
	    return jarPath.replaceAll(basePathName + ".*" ,basePathName + "/" + resourse);
	  }
	  
	  public static String getOs()
	  {
	    String os = System.getProperty("os.name").toLowerCase();
	    if (os.indexOf("win") >= 0)
	    {
	      return ("windows");
	    }
	    else if (os.indexOf("mac") >= 0)
	    {
	      return ("mac");
	    }
	    else if (os.indexOf("nix") >= 0 || os.indexOf("nux") >= 0 || os.indexOf("aix") > 0)
	    {
	      return ("linux");
	    }
	    else
	    {
	      throw new RuntimeException("Could not determine operating system");
	    }
	  }
	public static String getLine(String fileName, String searchString) 
	{
		try {
			BufferedReader reader = FileUtils.createReader(searchString);
			String line = null;
			while((line=reader.readLine())!=null)
				if(line.contains(searchString))
					return line;
			reader.close();
		}catch (FileNotFoundException e) {e.printStackTrace();} catch (IOException e) {e.printStackTrace();}
		return null;
	}
	public static String removeExtention(String jsonFN)
	{
		String[] eles = jsonFN.split("\\."); 
		String extention = "."+eles[eles.length-1];
		String withoutExtention = jsonFN.replace(extention, "");
		return withoutExtention;
	}
	
	public static String addBeforeExtention(String jsonFN, String addString) 
	{
		String newString = null;
		if(jsonFN.contains("."))
		{
			String[] eles = jsonFN.split("\\.");
			String extention = "."+eles[eles.length-1];
			String withoutExtention = jsonFN.replace(extention, "");
			newString = withoutExtention+addString+extention;
		}
		else
		{
			System.out.println("Warning, extension is missing on file:\t" + jsonFN);
			newString=jsonFN+addString;
		}
		return newString;
	}
	public static String makeFolderNameEndWithSlash(String fn) 
	{
		if(fn == null)
			return null;
		if(!fn.endsWith("/") && !fn.endsWith("\\"))
			fn = fn+"/";
		
		return fn;
	}
}
