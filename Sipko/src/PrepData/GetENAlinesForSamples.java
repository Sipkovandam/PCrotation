package PrepData;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;

public class GetENAlinesForSamples 
{
	//gets the ENA lines from an ENA file for those samples defined in samplesFN (should not contain anything but those names separated by enters).
	public static void main (String[] args) throws IOException
	{
		String ENA_FN = "E:/Groningen/Data/Juha/Calculon/JuhaMerged/2016-07-12-Q2-2016-ENA.txt";
		String samplesFN = "E:/Groningen/Data/Juha/Calculon/JuhaMerged/test/Missing.txt";
		String writeFN = "E:/Groningen/Data/Juha/Calculon/KallistoOnly2/2016-07-12-Q2-2016-ENA_DownloadsNeeded.txt";
		
		HashSet<String> samples = new HashSet<String>();
		BufferedReader reader = new BufferedReader(new FileReader(new File(samplesFN)));
		
		String line = null;
		
		while((line=reader.readLine())!=null)
			samples.add(line);
		reader.close();
		
		reader = new BufferedReader(new FileReader(new File(ENA_FN)));
		line = reader.readLine();//header line
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(writeFN)));
		System.out.println(line);
		writer.write(line+"\n");
		while((line=reader.readLine())!=null)
		{
			String[] cells = line.split("\t");
			String sample = cells[5];
			//System.out.println("sample =" + sample);
			if(samples.contains(sample))
				writer.write(line+"\n");
		}
		reader.close();
		writer.close();
		System.out.println("Done! File written to: " + writeFN);
	}
}
