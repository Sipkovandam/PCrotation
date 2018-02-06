package Tests;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.util.HashMap;
import java.util.concurrent.atomic.AtomicInteger;

import Slurm.ClusterVariables;
import Slurm.SharkVariables;
import Slurm.SlurmVariables;
import Tools.Script;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class Test extends Script<Test> {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public void test(double value)
	{
		HashMap<Double,Double> th = new HashMap<Double,Double>();
		th.put(value, value);
		System.out.println(value);
	}
	
	public void run()
	{
		
		String equal = "test";
		ClusterVariables slurmVariables = new SharkVariables();
		System.out.println("slurmVariables=" + slurmVariables.getMaxMemory());

		String test = "/groups/umcg-gdio/tmp04/umcg-svandam/Data/RNAseq/Counts/GRCh37/BloodStudy_BiosQcAnnique/Kallisto/Results/161221_NB501043_0096_AH775HBGX2_L1_ATGAGC"; 
		
		if(equal.equals("test"))
			System.out.println("test result:" + test.matches("_0694_|_0713_|_0096_"));
		
/*		final SamReader reader = SamReaderFactory.makeDefault().open(new File("bla"));
		for (final SAMRecord samRecord : reader) {
            // Convert read name to upper case.
            samRecord.setReadName(samRecord.getReadName().toUpperCase());
            samRecord.getCigar();
//            outputSam.addAlignment(samRecord);
        }
		
		long test = Long.parseLong("112345678910");
		long test2=0;
		test2 += (long)test;
		test2 += (long)test;
		p(test2);
		
		DecimalFormat f = new DecimalFormat("0000000000.0E0");

		double x = 123;
		if(x == Math.rint(x))
			p("bla");
		double val =1000000.06345345346;
		System.out.println(String.format("%.3f",val));
		System.out.println(String.format("%s = %d", "joe", 35));
		
		try
		{
			String gtfFn = "E:/Groningen/Data/Annotation/GRCh37/Homo_sapiens.GRCh37.75.gtf";
			GtfReader gtfReader = new GtfReader(new File(gtfFn));
			PerChrIntervalTree<GffElement> genome = gtfReader.createIntervalTree();
			List<GffElement> gffElements = genome.searchPosition("1", 11869);
			for(GffElement gffElement: gffElements)
			{
				System.out.println(gffElement.toString());
			}
			
		}catch(Exception e){e.printStackTrace();}
		
		
		PizzlyFusionStructure pizzlyFusionFile = new PizzlyFusionStructure();
		PizzlyFusionStructure pizzlyFusionStruct = (PizzlyFusionStructure) pizzlyFusionFile.readVars("E:/Groningen/Test/PizzlyClasses.PizzlyFusionFile/pizzly.json");
		p(pizzlyFusionStruct);
		System.out.println(pizzlyFusionStruct.getFusions());
		System.out.println(pizzlyFusionStruct.getFusions()[0].getGeneA());
		
		
		String n = "1	1635677	rs186584733	G	A	421.32	PASS	ABHom=0.998;AC=2;AF=1;AN=2;ANN=A|synonymous_variant|LOW|CDK11A|ENSG00000008128|transcript|ENST00000378633|protein_coding|15/20|c.1680C>T|p.His560His|1760/2458|1680/2352|560/783||,A|non_coding_transcript_exon_variant|MODIFIER|RP1-283E3.8|ENSG00000268575|transcript|ENST00000598846|processed_transcript|16/20|n.4289C>T||||||;BaseCounts=23,0,593,0;BaseQRankSum=-1.114;CADD=1.85435;CADD_SCALED=15.31;DB;DP=18;EXAC_AC_HET=148;EXAC_AC_HOM=6;EXAC_AF=0.001332;ExcessHet=0.0485;FS=0;GoNL_AF=0.0030120481927710845;GoNL_GTC=495|3|0;InbreedingCoeff=0.9999;MQ=53.36;MQRankSum=3.619;OND=0.002288;QD=19.15;ReadPosRankSum=2.032;SOR=0.914;VariantType=SNP;set=variant;RLV_PRESENT=TRUE;RLV=A|0.001332|CDK11A|0.0 0.02476038338658147|ENST00000378633||||||20164002_768974_DNA093895_445436_5GPM1610:HOMOZYGOUS||20164002_768974_DNA093895_445436_5GPM1610:1s1||Predicted pathogenic|GAVIN|Variant MAF of 0.001332 is rare enough to be potentially pathogenic and its CADD score of 15.31 is greater than a global threshold of 15.||;RLV_ALLELE=[A|CDK11A]A;RLV_ALLELEFREQ=[A|CDK11A]0.001332;RLV_GENE=[A|CDK11A]CDK11A;RLV_FDR=[A|CDK11A]0.0 0.02476038338658147;RLV_TRANSCRIPT=[A|CDK11A]ENST00000378633;RLV_PHENOTYPE=[A|CDK11A]NA;RLV_PHENOTYPEINHERITANCE=[A|CDK11A]NA;RLV_PHENOTYPEONSET=[A|CDK11A]NA;RLV_PHENOTYPEDETAILS=[A|CDK11A]NA;RLV_PHENOTYPEGROUP=[A|CDK11A]NA;RLV_SAMPLESTATUS=[A|CDK11A]20164002_768974_DNA093895_445436_5GPM1610:HOMOZYGOUS;RLV_SAMPLEPHENOTYPE=[A|CDK11A]NA;RLV_SAMPLEGENOTYPE=[A|CDK11A]20164002_768974_DNA093895_445436_5GPM1610:1s1;RLV_SAMPLEGROUP=[A|CDK11A]NA;RLV_VARIANTSIGNIFICANCE=[A|CDK11A]Predicted pathogenic;RLV_VARIANTSIGNIFICANCESOURCE=[A|CDK11A]GAVIN;RLV_VARIANTSIGNIFICANCEJUSTIFICATION=[A|CDK11A]Variant MAF of 0.001332 is rare enough to be potentially pathogenic and its CADD score of 15.31 is greater than a global threshold of 15.;RLV_VARIANTCOMPOUNDHET=[A|CDK11A]NA;RLV_VARIANTGROUP=[A|CDK11A]NA	GT:AD:DP:GQ:MQ0:PL	1/1:0,22:18:53:0:475,53,0";
		String split = n.split("(?=\\|ENSG)")[2].split("\\|")[1];
		//System.out.println(Arrays.toString(split));
		System.out.println(split);
		
		MyMatrix matrix1 = new MyMatrix("E:/Groningen/Test/MergeMatrixes/Matrix1.txt");
		MyMatrix matrix2 = new MyMatrix("E:/Groningen/Test/MergeMatrixes/Matrix2.txt");
		MyMatrix matrix3 = matrix1.mergeColumns(matrix2,2);
		matrix1.print();
		p("");
		matrix2.print();
		p("");
		matrix3.print();

		
		String test = "1";
		if(test.equals("1"))
			p("bla1");
		
		File file = new File("/Volumes/Promise_RAID/sipko/Data/Patrick/Merged/pcaFromRaw/counts_4309FineSamples_GPMbloedMerged_ScaffoldGenesRemoved_DiscardedRemoved_RawCounts.txt");
		String t = "AD10W1ACXX-7-23_1837_internal_id";
		String removed = FileUtils.removeExtention(t);
		System.out.println(removed);;
				
		String test = "bla";
		String test2 = "bla";
		
		if(test.equals(test2))
			System.out.println("ja");
		
		TestObject3<TestObject> test = new TestObject3<>();
		try {
			Thread.sleep(4123);
		} catch (InterruptedException e) {e.printStackTrace();}
		p("runtime = \t" + test.getRunTime());
		p(test.toString());
		*/
	}
	public void changeValue(int[] testVal)
	{
		testVal[0]= 10;
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
		a.log(a.jsonFN);
		
		return a;
	}

}

