package Tools;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
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
import java.io.RandomAccessFile;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.net.URLDecoder;

import org.apache.commons.io.IOUtils;

import PCA.RlogLargeMatrix;
import RowAnalyses.Row;

public class FileUtils
{
	public static HashMap<String, Double> readDoublehash(String fileName) throws FileNotFoundException, IOException
	{
		return readDoublehash(	fileName,
								0,
								1);
	}
	public static HashMap<String, Double> readDoublehash(String fileName, boolean hasHeader) throws FileNotFoundException, IOException
	{
		return readDoublehash(	fileName,
								0,
								1, hasHeader);
	}

	public static HashMap<String, Double> readDoublehash(	String fileName,
															int col1,
															int col2) throws FileNotFoundException, IOException
	{
		return readDoublehash(fileName,
												col1,
												col2, false);
	}
	
	public static void extractFile(String jarFn, String fileToGet)
	{
		String command = "jar xf " + jarFn + " " + fileToGet;
		System.out.println("Command =\t" + command);
		ExecCommand exec = new ExecCommand(command);
	}
	
	public static HashMap<String, Double> readDoublehash(	String fileName,
															int col1,
															int col2, boolean hasHeader) throws FileNotFoundException, IOException
	{
		HashMap<String, Double> hash = new HashMap<String, Double>();
		BufferedReader reader = createReader(fileName);
		String line = null;
		
		if(hasHeader)
			reader.readLine();
		
		while ((line = reader.readLine()) != null)
		{
			String[] cells = line.split("\t");
			hash.put(	cells[col1],
						Double.parseDouble(cells[col2].replace(	",",
																"")));
		}
		return hash;
	}

	public static ArrayList<String> readArrayList(String fileName) throws FileNotFoundException, IOException
	{
		return readArrayList(	fileName,
								false);
	}

	public static ArrayList<String> readArrayList(	String fileName,
													boolean firstColumnOnly) throws FileNotFoundException, IOException
	{
		if (fileName == null || !new File(fileName).exists())
		{
			throw new FileNotFoundException("File does not exist: " + fileName);
		}
		ArrayList<String> lines = new ArrayList<String>();

		BufferedReader reader = null;

		if (fileName.endsWith(".gz"))
			reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fileName))));
		else
			reader = new BufferedReader(new FileReader(new File(fileName)));

		String line = null;
		while ((line = reader.readLine()) != null)
		{
			if (firstColumnOnly)
				line = line.split("\t")[0];
			lines.add(line);
		}
		return lines;
	}
	public static BufferedReader createReader(String fn) throws FileNotFoundException, IOException
	{
		return createReader(fn, 8192);
	}
	

	public static BufferedReader createReader(String fn, int bufferSize) throws FileNotFoundException, IOException
	{
		BufferedReader reader = null;
		if (fn.endsWith(".gz"))
			reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fn))), bufferSize);
		else
			reader = new BufferedReader(new FileReader(new File(fn)), bufferSize);
		return reader;
	}

	
	public static BufferedWriter createWriter(String shellFN) throws FileNotFoundException, IOException
	{
		return createWriter(shellFN,
							false);
	}

	public static BufferedWriter createWriter(	String shellFN,
												boolean append) throws FileNotFoundException, IOException
	{
		BufferedWriter writer = null;
		if (shellFN.endsWith(".gz"))
			writer = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(	shellFN,
																											append))));
		else
			writer = new BufferedWriter(new FileWriter(	new File(shellFN),
														append));
		return writer;
	}

	public static void makeDir(String folderName)
	{
		if (folderName == null)
			return;

		File folder = new File(folderName);
		if (!folder.exists())
		{
			System.out.println("Folder does not exist, creating directory: \n" + folder);
			folder.mkdirs();
		}
	}

	public static String replaceEnd(String string,
									String addition)
	{
		String result = removeExtention(string) + addition;
		return result;
	}

	public static void copy(String kallistoShellFN,
							String tempKallistoFN) throws FileNotFoundException, IOException
	{
		ClassLoader.getSystemClassLoader();

		ClassLoader loader = FileUtils.class.getClassLoader();
		System.out.println(loader.getResource("foo/FileUtils.class"));

		InputStream inputStream = ClassLoader.getSystemResourceAsStream(kallistoShellFN);
		System.out.println("?=" + new File(kallistoShellFN).exists());
		System.out.println(inputStream);
		System.out.println(kallistoShellFN);
		InputStreamReader streamReader = new InputStreamReader(	inputStream,
																"UTF-8");
		BufferedReader reader = new BufferedReader(streamReader);
		BufferedWriter writer = createWriter(tempKallistoFN);
		String line = null;
		while ((line = reader.readLine()) != null)
		{
			System.out.println(line);
			writer.write(line + "\n");
		}
	}
	public static HashMap<String, Integer> makeHash(String line)
	{
		return makeHash(line, "\t");
	}

	public static HashMap<String, Integer> makeHash(String line, String sep)
	{
		String[] cells = line.split(sep);
		HashMap<String, Integer> index = new HashMap<String, Integer>();
		for (int c = 0; c < cells.length; c++)
		{
			index.put(	cells[c],
						c);
		}
		return index;
	}

	public static HashMap<String, String> readStringStringHash(String ensgToGeneSymbolFN)
	{
		return readStringStringHash(ensgToGeneSymbolFN,
									false,
									0);
	}
	
	public static HashMap<String, String> readStringStringHash(String ensgToGeneSymbolFN, int skipLines)
	{
		return readStringStringHash(ensgToGeneSymbolFN,
									false,
									skipLines);
	}

	public static HashMap<String, String> readStringStringHash(	String ensgToGeneSymbolFN,
																int keyCol,
																int valueCol)
	{
		return readStringStringHash(ensgToGeneSymbolFN,
									false,
									0,
									keyCol,
									valueCol);
	}

	public static HashMap<String, String> readStringStringHash(	String ensgToGeneSymbolFN,
																boolean wholeLineValue,
																int skipHeaderRows)
	{
		return readStringStringHash(ensgToGeneSymbolFN,
									wholeLineValue,
									skipHeaderRows,
									0,
									1);
	}

	public static HashMap<String, String> readStringStringHash(	String ensgToGeneSymbolFN,
																boolean wholeLineValue,
																int skipHeaderRows,
																int keyCol,
																int valueCol)
	{
		if(!new File(ensgToGeneSymbolFN).exists())
		{
			System.out.println("Warning file does not exist:\t" + ensgToGeneSymbolFN);
			return null;
		}

		
		HashMap<String, String> conversionHash = new HashMap<String, String>();
		try
		{
			BufferedReader reader = createReader(ensgToGeneSymbolFN);
			reader.lines().skip(skipHeaderRows).forEach(line ->
			{
				String[] cells = line.split("\t",
											2);
				if (cells.length > 1 && wholeLineValue)
					conversionHash.put(	cells[keyCol],
										cells[1]);
				else if (!wholeLineValue)
				{
					cells = line.split("\t");
					if(cells.length==1)
						conversionHash.put(	cells[keyCol],
						                   	"");
					else
						conversionHash.put(	cells[keyCol],
					                   	cells[valueCol]);
				}
			});
		} catch (FileNotFoundException e)
		{
			e.printStackTrace();
		} catch (IOException e)
		{
			e.printStackTrace();
		}
		return conversionHash;
	}

	public static HashMap<String, Double> makeHash(	String readCountFile,
													int i)
	{
		HashMap<String, Double> keyToValue = new HashMap<String, Double>();
		try
		{
			BufferedReader reader = createReader(readCountFile);
			reader.lines().forEach(line ->
			{
				String[] eles = line.split("\t");
				keyToValue.put(	eles[0],
								Double.parseDouble(eles[i]));
			});
			reader.close();
		} catch (IOException e)
		{
			e.printStackTrace();
		}
		return keyToValue;
	}

	public static InputStream getResourceFile(String file) throws IOException
	{
		return Util.class.getResourceAsStream("/" + file);
	}

	public static String getResourcePath(String file) throws IOException
	{
		URL url = Util.class.getResource("/" + file);
		if (url == null)
			throw new IOException("resource file not found: " + file);
		return url.toString().replaceAll(	"file:/",
											"file:///");
	}

	public static String readFile(String file) throws IOException
	{
		return IOUtils.toString(Util.getResourceFile(file));
	}

	public static String prepareFileFromJar(String file) throws IOException
	{
		return FileUtils.prepareFileIntFromJar(	file,
												false,
												true);
	}

	public static String prepareBinaryFromJar(String file) throws IOException
	{
		return FileUtils.prepareFileIntFromJar(	file,
												true,
												true);
	}

	public static String getTempFile(String file) throws IOException
	{
		return FileUtils.prepareFileIntFromJar(	file,
												false,
												false);
	}

	public static void cleanTemp() throws IOException
	{
		String tempPath = FileUtils.getTempPath();
		(new File(tempPath)).delete();
	}

	private static String prepareFileIntFromJar(String file,
												boolean isBinary,
												boolean doCreate) throws IOException
	{
		// Copies to location next to JAR file
		String tempPath = FileUtils.getTempPath();
		(new File(tempPath)).mkdir();
		File outBinFileName = new File(tempPath + file.substring(Math.max(	0,
																			file.lastIndexOf("/"))));
		if (doCreate)
		{
			String fn = file.replaceAll("\\{OS\\}",
										getOs());
			System.out.println(fn + "\t" + new File(fn).exists() + "\t" + outBinFileName);

			InputStream stream = Util.getResourceFile(fn);
			FileOutputStream output = new FileOutputStream(outBinFileName);
			IOUtils.copy(	stream,
							output);
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
		String jarPath = URLDecoder.decode(	path,
											"UTF-8").toString();
		return jarPath.substring(	0,
									jarPath.lastIndexOf("/"))
				+ "/temp/";
	}

	public static String getBigResource(Class clazz,
										String basePathName,
										String resourse) throws IOException
	{
		String path = clazz.getProtectionDomain().getCodeSource().getLocation().getPath();
		String jarPath = URLDecoder.decode(	path,
											"UTF-8").toString();
		return jarPath.replaceAll(	basePathName + ".*",
									basePathName + "/" + resourse);
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

	public static String getLine(	String fileName,
									String searchString)
	{
		try
		{
			BufferedReader reader = FileUtils.createReader(fileName);
			String line = null;
			while ((line = reader.readLine()) != null)
				if (line.contains(searchString))
					return line;
			reader.close();
		} catch (FileNotFoundException e)
		{
			e.printStackTrace();
		} catch (IOException e)
		{
			e.printStackTrace();
		}
		return null;
	}

	public static String removeExtention(String fileName)
	{
		fileName = fileName.replace(".gz",
									"");//if it is zipped remove this bit first
		String[] eles = fileName.split("\\.");
		String extention = "." + eles[eles.length - 1];
		String withoutExtention = fileName.replace(	extention,
													"");
		return withoutExtention;
	}
	
	public static double[] convertToDoubleArray(String[] valuesString)
	{
		return convertToDoubleArray(valuesString, new double[valuesString.length]);
	}
	
	public static double[] convertToDoubleArray(String[] valuesString, double[] values)
	{
		for (int v = 0; v < valuesString.length; v++)
			values[v] = Double.parseDouble(valuesString[v]);
		return values;
	}

	public static String addBeforeExtention(String fileName,
											String addString)
	{
		String newString = null;
		if (fileName.contains("."))
		{
			boolean isZipped = false;
			if (fileName.endsWith(".gz"))
			{
				String[] eles = fileName.split("\\.");
				String extention = "." + eles[eles.length - 1];
				String withoutExtention = fileName.replace(	extention,
															"");
				fileName = withoutExtention;
				isZipped = true;
			}

			String[] eles = fileName.split("\\.");
			String extention = "." + eles[eles.length - 1];
			String withoutExtention = fileName.replace(	extention,
														"");
			newString = withoutExtention + addString + extention;

			if (isZipped)
			{
				newString += ".gz";
			}
		}
		else
		{
			System.out.println("Warning, extension is missing on file:\t" + fileName);
			newString = fileName + addString;
		}
		return newString;
	}

	public static String makeFolderNameEndWithSlash(String fn)
	{
		if (fn == null)
			return null;
		if (!fn.endsWith("/") && !fn.endsWith("\\"))
			fn = fn + "/";

		return fn;
	}

	public static HashMap<String, Integer> StringToIndexHash(String fn1)
	{
		return StringToIndexHash(	fn1,
									true);
	}

	public static HashMap<String, Integer> StringToIndexHash(	String fn1,
																boolean header)
	{
		HashMap<String, Integer> rowNameToIndex = new HashMap<String, Integer>();
		try
		{
			BufferedReader reader = createReader(fn1);
			if (header)
				reader.readLine();
			String line = null;
			int i = 0;
			while ((line = reader.readLine()) != null)
			{
				String rowName = line.split("\t")[0];
				if (!rowNameToIndex.containsKey(rowName))
					rowNameToIndex.put(	rowName,
										i);
				else
				{
					System.out.println("Duplicate rowNames:" + rowName);
					System.out.println("Exiting");
					System.exit(2);
				}
				i++;
			}

		} catch (FileNotFoundException e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return rowNameToIndex;
	}

	public static HashMap<String, ArrayList<String>> readStringMultiStringArrayList(	String conversionFile,
																				int i,
																				int j,
																				boolean wholeLineValue)
	{
		HashMap<String, ArrayList<String>> conversionHash = new HashMap<String, ArrayList<String>>();
		try
		{
			BufferedReader reader = createReader(conversionFile);
			String line = null;
			while((line=reader.readLine())!=null)
			{
				String[] cells = line.split("\t",2);
				if (cells.length > 1 && wholeLineValue)
				{
					conversionHash=addToMultiKeyHash(cells[i], cells[1], conversionHash);

				}
				else if (cells.length > 1 && !wholeLineValue)
				{
					cells=line.split("\t");
					conversionHash=addToMultiKeyHash(cells[i], cells[j], conversionHash);
				}
			}
		} catch (FileNotFoundException e)
		{
			e.printStackTrace();
		} catch (IOException e)
		{
			e.printStackTrace();
		}
		return conversionHash;
	}
	public static HashMap<String, HashMap<String, Integer>> readStringMultiStringHash(	String conversionFile,
	    																				int i,
	    																				int j,
	    																				boolean wholeLineValue)
	{
		return readStringMultiStringHash(	conversionFile,
											i,
											j,
											wholeLineValue, false);
	}
	
	public static HashMap<String, HashMap<String, Integer>> readStringMultiStringHash(	String conversionFile,
																				int i,
																				int j,
																				boolean wholeLineValue, boolean skipFirstline)
	{
		HashMap<String, HashMap<String, Integer>>  conversionHash = new HashMap<String, HashMap<String, Integer>>();
		try
		{
			BufferedReader reader = createReader(conversionFile);
			if(skipFirstline)
				reader.readLine();
			String line = null;
			while((line=reader.readLine())!=null)
			{
				String[] cells = line.split("\t",2);
				if (cells.length > 1 && wholeLineValue)
				{
					conversionHash=addToMultiKeyArrayList(cells[i], cells[1], conversionHash);

				}
				else if (cells.length > 1 && !wholeLineValue)
				{
					cells=line.split("\t");
					conversionHash=addToMultiKeyArrayList(cells[i], cells[j], conversionHash);
				}
			}
		} catch (FileNotFoundException e)
		{
			e.printStackTrace();
		} catch (IOException e)
		{
			e.printStackTrace();
		}
		return conversionHash;
	}

	private static  HashMap<String, HashMap<String, Integer>> addToMultiKeyArrayList(String key,
																		String value,
																		HashMap<String, HashMap<String, Integer>> conversionHash)
	{
		HashMap<String, Integer> values = conversionHash.get(key);
		if(values == null)
			values = new HashMap<String, Integer>();
		values.put(value,0);
		conversionHash.put(key, values);
		return conversionHash;
	}
	private static HashMap<String, ArrayList<String>> addToMultiKeyHash(String key,
																		String value,
																		HashMap<String, ArrayList<String>> conversionHash)
	{
		ArrayList<String> values = conversionHash.get(key);
		if(values == null)
			values = new ArrayList<String>();
		values.add(value);
		conversionHash.put(key, values);
		return conversionHash;
	}

	public static Row readRow(String line)
	{
		return Row.readRow(line);
	}

	public static HashMap<String, ArrayList<String>> readGavinVcfHash(String vcfFn, String enstToEnsgFn) throws FileNotFoundException, IOException
	{
		if(!new File(vcfFn).exists())
		{
			System.out.println("Warning file does not exist:\t" + vcfFn);
			return null;
		}
		
		HashMap<String,String> enstToEnsg = readStringStringHash(enstToEnsgFn);
		
		BufferedReader vcfReader = createReader(vcfFn);
		HashMap<String, ArrayList<String>> ensgToPositions = new HashMap<String, ArrayList<String>>();
		String line = null;
		while((line=vcfReader.readLine())!=null)
		{
			String transcriptId = line.split("\t")[10];
			String ensgId = enstToEnsg.get(transcriptId);
			if(transcriptId.contains("ENSG"))
				ensgId=transcriptId;
			if(ensgId == null)
			{
//				if(!transcriptId.contains("Transcript"))//header
//					System.out.println("Warning no EnsemblId could be found for transcript: " + transcriptId);
				continue;
			}

			ArrayList<String> positions = ensgToPositions.get(ensgId);
			if(positions == null)
				positions=new ArrayList<String>();

			positions.add(line);
			ensgToPositions.put(ensgId, positions);
		}
		return ensgToPositions;
	}

	public static String getMd5(String countFilePath) throws IOException
	{
		FileInputStream fis = new FileInputStream(new File(countFilePath));
		String md5 = org.apache.commons.codec.digest.DigestUtils.md5Hex(fis);
		fis.close();
		return md5;
	}
	
	public static void checkNull(Object pizzlyExecute, String string)
	{
		if(pizzlyExecute == null)
		{
			System.out.println(string+" variable is not initiated; Exiting");
			System.exit(2);
		}
	}

	public static void writeHashtable(HashMap<String, int[]> countsPerString, String writeFn, String header) throws FileNotFoundException, IOException
	{
		BufferedWriter writer = createWriter(writeFn);
		
		if(header != null)
			writer.write(header+"\n");
		
		for(String name : countsPerString.keySet())
		{
			String line = name;
			int[] values = countsPerString.get(name);
			for(int value : values)
			{
				line = line.concat("\t").concat(Integer.toString(value));
			}
			line=line.concat("\n");
			writer.write(line);
		}
		
		writer.close();
	}
	public static HashMap<String, Integer> arrayToHashMap(String[] strings)
	{
		HashMap<String, Integer> stringToIndex = new HashMap<String, Integer>();
		for(int s = 0; s < strings.length; s++)
		{
			stringToIndex.put(strings[s], s);
		}
		return stringToIndex;
	}
	public static String doubleArrayToWriteString(double[] values)
	{
		String line = "";
		StringBuilder stringbuilder = new StringBuilder();
		for(double value : values)
		{
			stringbuilder.append("\t").append(value);
		}
		line=stringbuilder.toString();
		return line;
	}
	public static String StringArrayToWriteString(String[] dataColHeaders)
	{
		String line = "";
		StringBuilder stringbuilder = new StringBuilder();
		for(String dataColHeader : dataColHeaders)
		{
			stringbuilder.append("\t").append(dataColHeader);
		}
		line=stringbuilder.toString();
		return line;
	}
	public static HashMap<String, Integer> readStringToIntegerHash(String indexFile) throws NumberFormatException, IOException
	{
		HashMap<String, Integer> hash = new HashMap<String, Integer>();
		BufferedReader reader = createReader(indexFile);
		String line = null;
		
//		if(indexFile)
//			reader.readLine();
		
		while ((line = reader.readLine()) != null)
		{
			String[] cells = line.split("\t");
			hash.put(	cells[0],
			         	Integer.parseInt(cells[1].replace(	",",
																"")));
		}
		return hash;
	}
	public static HashMap<String, Long> readStringToLongHash(String indexFn) throws FileNotFoundException, IOException
	{
		HashMap<String, Long> hash = new HashMap<String, Long>();
		BufferedReader reader = createReader(indexFn);
		String line = null;
		
//		if(indexFile)
//			reader.readLine();
		
		while ((line = reader.readLine()) != null)
		{
			String[] cells = line.split("\t");
			
			if(hash.containsKey(cells[0]))
			{
				System.out.println("File contains the same rowName multiple times, recommended to use readStringMultiLongArrayList instead");
			}
			
			hash.put(	cells[0],
			         	Long.parseLong(cells[1].replace(	",",
																"")));
		}
		return hash;
	}
	public static String getFolderName(String fn)
	{
		return new File(new File(fn).getParent()).getName();
	}
	public static String getFolderName(File fn)
	{
		return new File(fn.getParent()).getName();
	}
	public static ArrayList<String> addStringToArrayList(ArrayList<String> names, String name)
	{
		if(names==null)
			names = new ArrayList<String>();
		names.add(name);
		return names;
	}
	
	public static BufferedReader getWebsiteReader(String webUrl) throws IOException
	{
		URL geneNetwork = new URL(webUrl);
		System.out.println(webUrl);
        URLConnection yc =geneNetwork.openConnection();
        System.out.println("yc=" + yc);
    	BufferedReader in = null;
        try
        {
        	in = new BufferedReader(new InputStreamReader(yc.getInputStream(), "UTF-8"));
        }catch(Exception e)
        {
        	e.printStackTrace();
        	System.out.println("Errro cannot open this web url: \t" + webUrl);
        	System.exit(1);
        }
        return in;
	}
	public static void saveFileFromUrl(String url, String fn) throws IOException
	{
		BufferedInputStream in = new BufferedInputStream(new URL(url).openStream());
		BufferedOutputStream writer = new BufferedOutputStream(new FileOutputStream(fn));
	    byte data[] = new byte[1024];
	    int count;
	    while((count = in.read(data,0,1024)) != -1)
	    {
	    	writer.write(data, 0, count);
	    }
	    writer.close();
	}
	public static String makeLineFromArray(String[] colValues)
	{
		if(colValues.length==0)
			return null;
		
		StringBuilder lineBuilder = new StringBuilder();
		lineBuilder.append(colValues[0]);
		for(int c =1; c <  colValues.length; c++)
		{
			lineBuilder.append("\t");
			String toAdd =colValues[c];
			if(toAdd==null)
				toAdd="";
			lineBuilder.append(toAdd);
		}
		
		return lineBuilder.toString();
	}
}
