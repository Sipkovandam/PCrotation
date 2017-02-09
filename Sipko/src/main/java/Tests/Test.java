package Tests;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
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
import java.lang.reflect.*;

import org.junit.rules.TemporaryFolder;

import PCA.Matrix;
import PCA.MatrixStruct;
import PCA.RLog;
import PCA.RlogLargeMatrix;
import PCA_testcases.CompareFiles;
import STAR.SpliceSitesCount;
import STAR._STAR_Pipeline;
import Tools.FileUtils;
import Tools.Script;
import umcg.genetica.containers.Pair;
import umcg.genetica.math.matrix.DoubleMatrixDataset;


public class Test extends Script<Test> {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public void run()
	{
		String test = "bla";
		String test2 = "bla";
		
		if(test.equals(test2))
			System.out.println("ja");
		/*
		TestObject3<TestObject> test = new TestObject3<>();
		try {
			Thread.sleep(4123);
		} catch (InterruptedException e) {e.printStackTrace();}
		p("runtime = \t" + test.getRunTime());
		p(test.toString());
		*/
	}
	
	public static void main(String[] args)  throws IOException, CloneNotSupportedException, ClassNotFoundException, InstantiationException, IllegalAccessException, NoSuchMethodException, SecurityException, IllegalArgumentException, InvocationTargetException, InterruptedException
	{	
		new Test().run();
		
/*		String test = "a,123,asdf,n,niet,asfdasdf,1234,0,1,f,appel,kaas,123,1234\t999";
		String[] split =  test.split("(?<=\\G.*,.*,.*,.*,.*,.*),");
		Stream.of(split).forEach(Script::p);
		p(" ");

		String[] split2 = test.split("(?<=\\G.*,.*,.*,.*,.*),"); // 
		
		Stream.of(split2).forEach(Script::p);
		
		//TestObject3<TestObject> testObject3_2 = 

		String param1= "appel";
		Double param2= new Double(12);
		String className = "Tests.TestObject2";
		Class<?> cl = (Class<?>) Class.forName(className);
		Constructor con = cl.getDeclaredConstructor(String.class, Double.class);
		Object xyz = (Object) con.newInstance(param1, param2);
		
//		p(xyz.getClass().getName());
//		p(xyz.getClassName());
//		p(xyz.getExecuteClass());
		
//		TestObject2 x = (TestObject2) xyz;
//		
//		p(x.name);
//		
		for(int m = 0; m < cl.getMethods().length; m++)
		{
			p("" + cl.getMethods()[m].getName());
			if( cl.getMethods()[m].getAnnotatedParameterTypes().length>0)
			{
				p("" + cl.getMethods()[m].getAnnotatedParameterTypes()[0]);
				p("" + cl.getMethods()[m].getAnnotatedParameterTypes()[0]);
			}
		}
		
		Method method = cl.getMethod("setExecuteClass", String.class);
//
		method.invoke(xyz, "ja");
//
		Method method2 = cl.getMethod("getExecuteClass");
		p(""+method2.invoke(xyz));
		
		
		
		
//		String cn = "java.io.File";
//		Class myClass = Class.forName(cn);
//		
//		File testFile = new File("aap");
//		
//		Class[] types = {testFile.getClass()};
//		
//		Constructor constructor = myClass.getConstructor(types);
//		
//		Object[] parameters = {File.class};
//		Object instanceOfMyClass = constructor.newInstance(parameters);
//		
//		p(instanceOfMyClass.getClass().getName());
		
		
		/*String test = "a,123,asdf,n,niet,asfdasdf,1234,0,1,f,appel,kaas,123,1234\t999";
		String[] split =  test.split("(?<=\\G.*,.*,.*,.*,.*,.*),");
		Stream.of(split).forEach(Script::p);
		p(" ");
		String data = "a,123,asdf,n,niet,asfdasdf,1234,0,1,f,appel,kaas,123,1234,999,123";
		
		String[] split2 = data.split("(?<=\\G.*,.*,.*,.*,.*),"); //
		
		Stream.of(split2).forEach(Script::p);
		
		 Creating 100K readers is pretty fast. Yay!
		int nReaders=100000;
		BufferedReader[] reader = new BufferedReader[nReaders];
		
		for(int r = 0;r < nReaders; r++)
		{
			if(r%10000==0)
				p(r+"/"+nReaders);
			reader[r] = FileUtils.createReader("E:/Groningen/Test/STAR/STAR/Results/StartEndToGene.txt");
		}
		for(int r = 0;r < nReaders; r++)
		{
			System.out.println(reader[r].readLine());
			reader[r].close();
		}
		
		
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
		System.out.println("Test done!");
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
	
	private static Test t(Test a)
	{
		a.p(a.jsonFN);
		
		return a;
	}

}

