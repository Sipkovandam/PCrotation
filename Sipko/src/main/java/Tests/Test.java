package Tests;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import org.junit.rules.TemporaryFolder;

import PCA.Matrix;
import PCA.MatrixStruct;
import PCA.RLog;
import PCA.RlogLargeMatrix;
import PCA_testcases.CompareFiles;
import Tools.Script;
import umcg.genetica.containers.Pair;
import umcg.genetica.math.matrix.DoubleMatrixDataset;


public class Test {

	public static void main(String[] args) throws IOException, CloneNotSupportedException, ClassNotFoundException
	{	
		
		
		/*
		String s = "TestObject";

	    Object o = null;

	    System.out.println(o.getClass());
	    Class<? extends Script<RlogLargeMatrix>> type;
	    
	    try {
	        Class<?> type2 = Class.forName(s);
	        getValue(type2, o);
	    } catch (ClassNotFoundException e) {
	        // class was not found
	        e.printStackTrace();
	    }
		
		System.out.println(o);
		
		String testString = "appelflap";
		String testS2 = "taart";
		
		File file = new File("E:/Groningen/Splicing/100BPcap_analysis/STAR_config_200Samples_100BPcap_analysis.json");
		
		for(int i = 0; i<250000000;i++)
		{
			if(i%1000000==0)
				System.out.println("i=" + i +"/25000000");
			if(file.getName().contains(testString))
				System.out.println("bla");
		}
		
		
		
		double[] genes = new double[]{1,2,10,20};
		double geneCounts = Stream.of(genes).collect(Collectors.summarizingDouble(g -> 
		{
			System.out.println(g);
			return 0;
		})).getSum();
		
		
		ArrayList<String> sTAR_SpliceFiles = new ArrayList<String>();
		sTAR_SpliceFiles.add("asdf");
		sTAR_SpliceFiles.add("234");
		sTAR_SpliceFiles.add("fes2");
		sTAR_SpliceFiles.add("jkhl");
		String[] rowNames = new String[sTAR_SpliceFiles.size()];
		sTAR_SpliceFiles.toArray(rowNames);
		
		Stream.of(rowNames).forEach(System.out::println);
		
		
		TestObject var = new TestObject();
		var.bla="Neej!";
		var.ts = new String[]{"ja","toch"};
		System.out.println(var.bla);
		System.out.println(var.ts[0] + " " + var.ts[1]);
		var.hash = new HashMap<String,String>();
		var.hash.put("first", "firstVal");
		var.hash.put("2nd", "2ndVal");
		System.out.println(var.hash.get("first") + " " + var.hash.get("2nd"));
		
		TestObject copy = (TestObject) var.clone();
		System.out.println(copy.bla);
		System.out.println(copy.ts[0] + " " + copy.ts[1]);
		System.out.println(copy.hash.get("first") + " " + copy.hash.get("2nd"));
		
		HashMap<String, double[]> test = new HashMap<String, double[]>();
		test.put("bla",new double[2]);
		System.out.println(test.containsKey(null));
		
		
		Pair<Integer,Integer> test = new Pair<Integer,Integer>(3,5);
		
		String test = "bla\tbla\t-\t-\t\t\t\t";
		
		String[] cells = test.split("\t");
		
		System.out.println(cells.length);
		
		
		HashMap<String,String> test = new HashMap<String,String>();
		test.put("test", "out1");
		test.put("test2", "out2");
		test.put("test3", "out13");
		
		List<Entry<String, String>> testout = test.entrySet().stream().collect(Collectors.toList());
		testout.forEach(entry -> System.out.println(entry.getKey()));
		
		
		String test = "123";
		addToString(test);
		System.out.println(test);
		
		String folder = "/home/directory/test/";
		String[] folders = folder.split(",");
		System.out.println(folders[0]);
		
		File file = new File("/test/directory/parentget/");
		System.out.println(file.getParent());
		
		
		double testval= 1.7;
		
		double[] testVals = new double[]{1,2,3,4,5};
		IntStream.range(0, testVals.length).forEach(x -> System.out.println(testVals[x]));

		Hashtable<String,Integer> testhash = new Hashtable<String,Integer>();
		testhash.put("1", 1);
		testhash.put("2", 2);
		testhash.put("3", 3);
		testhash.put("4", 4);

		System.out.println("hashtable");
		testhash.forEach((key,value) -> System.out.println(key +" " + value));
		
		System.out.println("test");
		String line = "1	13053	13220	1	3	1	1	0	25";
		AtomicInteger tab = new AtomicInteger(0);
		AtomicInteger pos = new AtomicInteger(0);
		AtomicInteger splitChar = new AtomicInteger(0);
		line.chars().forEach(c -> {
			if(c== 9)//9 == tab
				tab.getAndIncrement();
				if(tab.get()==6)
					splitChar.set(pos.get());
				pos.getAndIncrement();
					});
		
		System.out.println("tab= " + splitChar.get() + " left =" + line.substring(0,splitChar.get()));
//		String geoFN= "E:/Groningen/Scripts/Tests/Rlog.java/DESeqNorm/Samples.DESeqNorm.txt.gz";
//		CompareFiles.compare(geoFN,"TestCaseFiles/Result_Samples.txt",true);
//		
*/
	}
	private static <T> T getValue(Class<T> desiredType, Object o) 
	{
		if (o.getClass().isAssignableFrom(desiredType)) 
		{
			return desiredType.cast(o);
		} else 
		{
			throw new IllegalArgumentException();
		}
    }
    
	private static void addToString(String testString)
	{
		testString+="blblbla";
		System.out.println(testString);
	}
	private static void testFunc(int c, AtomicInteger tab) 
	{
			if(c== 9)//9 == tab
				tab.getAndIncrement();
				if(tab.get()==6)
					System.out.println("bla");
	}
	public void print()
	{
		
		
	}
}

