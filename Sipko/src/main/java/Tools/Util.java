package Tools;

//Author: Jesse van dam
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.net.URLDecoder;

import org.apache.commons.io.IOUtils;

public class Util
{
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
    return Util.prepareFileIntFromJar(file,false,true);
  }
  
  public static String prepareBinaryFromJar(String file) throws IOException
  {
    return Util.prepareFileIntFromJar(file,true,true);
  }
  
  public static String getTempFile(String file) throws IOException
  {
    return Util.prepareFileIntFromJar(file,false,false);
  }
  
  public static void cleanTemp() throws IOException
  {
    String tempPath = Util.getTempPath();
    (new File(tempPath)).delete();
  }
  
  private static String prepareFileIntFromJar(String file, boolean isBinary, boolean doCreate) throws IOException
  {
    // Copies to location next to JAR file
    String tempPath = Util.getTempPath();
    (new File(tempPath)).mkdir();
    File outBinFileName = new File(tempPath + file.substring(Math.max(0,file.lastIndexOf("/"))));
    System.out.println("FilePath =" + file);
    System.out.println("FilePathExists =" + new File(file).exists());
    if (!outBinFileName.exists() && doCreate)
    {
      System.out.println(file.replaceAll("\\{OS\\}",getOs()) + "\n" + outBinFileName);
      InputStream stream = Util.getResourceFile(file.replaceAll("\\{OS\\}",getOs()));
      FileOutputStream output = new FileOutputStream(outBinFileName);
      System.out.println("File = " + file.replaceAll("\\{OS\\}",getOs()));
      System.out.println("File exists = " + new File(file.replaceAll("\\{OS\\}",getOs())).exists());
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
  
}