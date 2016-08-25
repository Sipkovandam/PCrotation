package PrepData;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashSet;

public class KallistoRemoveCompleted {
	

	public static void main(String[] ags) throws Exception
	{
		String enaFN = "E:/Groningen/Data/Juha/Calculon/StartKallisto4/2016-07-12-Q2-2016-ENA.txt";
		String completedFN = "E:/Groningen/Data/Juha/Calculon/StartKallisto4/completed2.txt";
		String writeFN = "E:/Groningen/Data/Juha/Calculon/StartKallisto4/2016-07-12-Q2-2016-ENA_left2.txt";
		
		HashSet<String> completeHash = new HashSet<String>();
		BufferedReader reader = new BufferedReader(new FileReader(new File(completedFN)));
		String line = null;
		while((line=reader.readLine())!=null)
		{
			String[] eles = line.split("/");
			System.out.println(eles[eles.length-3]);
			completeHash.add(eles[eles.length-3]);
		}
		System.out.println("Nubmer of completed files = " + completeHash.size());
		
		reader = new BufferedReader(new FileReader(new File(enaFN)));
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(writeFN)));
		line = null;
		while((line=reader.readLine())!=null)
		{
			String[] eles = line.split("\t");
			String sample = eles[5];
			if(completeHash.contains(sample))
				continue;
			writer.write(line+"\n");
		}
		System.out.println("Done, file written to:" + writeFN);
		writer.close();
		reader.close();
	}
}
