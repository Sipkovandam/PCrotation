package PrepData;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;

public class FailedFilesKallisto 
{
	public static void main (String[] args) throws IOException
	{
		String missingDownloads =new String("E:/Groningen/Data/Juha/Calculon/KallistoOnly/missingDownloads2.txt");//obtained with:find /groups/umcg-wijmenga/tmp04/umcg-svandam/Juha/KallistoRun/logs_out -name '*.out' -exec grep "'n_not_downloaded', '1'" {} \; -exec echo {} \; | grep groups > missingDownloads.txt
		String failedKallistoFiles =new String("E:/Groningen/Data/Juha/Calculon/KallistoOnly/kallistoFail2.txt");//obtained with:find /groups/umcg-wijmenga/tmp04/umcg-svandam/KallistoResults/KallistoResults/ -name kallisto_return_code.txt -exec grep -v "1" {} \; -exec echo {} \; > kallistoFail.txt
		String ENA_Lines = new String("E:/Groningen/Data/Juha/Calculon/KallistoOnly/2016-04-20-ENA_ALL_Q2-2016_AllNoDups.txt");
		String mappingPerc = new String("E:/Groningen/Data/Juha/Calculon/KallistoOnly/mappingPerSample2.txt");//obtained with:find /groups/umcg-wijmenga/tmp04/umcg-svandam/KallistoResults/KallistoResults/  -name '*.err' -exec grep -H "pseudoaligned" {} \; | tr ' ' '\t' | cut -f1,3,5 | sed -e 's/\[quant\]//g' > mappingPerSample.txt
		String writeFN = "E:/Groningen/Data/Juha/Calculon/KallistoOnly/2016-04-20-ENA_ALL_Q2-2016_AllNoDups_MappingInfoAdded.txt";
		
		String line = null;
		BufferedReader reader = new BufferedReader(new FileReader(new File(ENA_Lines)));
		int n = 0;
		Hashtable<Integer,String> ENA_Line_Hash = new Hashtable<Integer,String>();
		while((line = reader.readLine())!=null)
		{
			ENA_Line_Hash.put(n, line);
			n++;
		}
		reader.close();
		
		//put failed downloads in hash
		reader = new BufferedReader(new FileReader(new File(missingDownloads)));
		Hashtable<String,String> missing = new Hashtable<String,String>();
		while((line = reader.readLine())!=null)
		{
			int lineNumber = Integer.parseInt(line.replace("/groups/umcg-wijmenga/tmp04/umcg-svandam/Juha/KallistoRun/logs_out/","").replace(".out", ""));
			String ENA_Line = ENA_Line_Hash.get(lineNumber+1);
			System.out.println(lineNumber+"\t"+ENA_Line);
			String sample = ENA_Line.split("\t")[5];
			System.out.println(sample);
			missing.put(sample, "fastq download missing");
		}
		reader.close();
		
		//put failed Kallisto runs in hash
		reader = new BufferedReader(new FileReader(new File(failedKallistoFiles)));
		Hashtable<String,Integer> failedKallisto = new Hashtable<String,Integer>();
		
		int lineNumber = 1;
		int failCode = 0;
		System.out.println("");
		while((line = reader.readLine())!=null)
		{
			if(lineNumber%2==1)
			{
				failCode = Integer.parseInt(line);
				lineNumber++;
				continue;
			}
			
			String sample = line.split("/")[9];
			failedKallisto.put(sample, failCode);
			lineNumber++;
			failCode = 0;
		}
		reader.close();
		
		//put aligned reads in 2 hashes
		reader = new BufferedReader(new FileReader(new File(mappingPerc)));
		Hashtable<String,double[]> mappingPerc29 = new Hashtable<String,double[]>();
		Hashtable<String,double[]> mappingPerc31 = new Hashtable<String,double[]>();
		
		while((line = reader.readLine())!=null)
		{			
			String sample = line.split("/")[9];
			///System.out.println(sample);
			double total= Double.parseDouble(line.split("\t")[1].replace(",", ""));
			double mapping= Double.parseDouble(line.split("\t")[2].replace(",", ""));
			double percentage= mapping/total;
			double[] values = new double[]{mapping,total,percentage};
			if(line.contains("/29/"))
			{
				mappingPerc29.put(sample, values);
			}
			if(line.contains("/31/"))
			{
				mappingPerc31.put(sample, values);
			}
		}
		reader.close();
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(writeFN)));
		reader = new BufferedReader(new FileReader(new File(ENA_Lines)));
		
		String headerLine ="Reads_Mapping_29\tReads_Mapping_31\tFold_Difference\tReads_total\tPercentage_mapping29\tPercentage_mapping31\tmissingDownload\tKallisto_return_code\t"+reader.readLine()+"\n";
		writer.write(headerLine);
		while((line = reader.readLine())!=null)
		{
			String sample = line.split("\t")[5];
			//System.out.println(sample);
			
			String outLine="";
			double mapping29 = Double.NaN;
			double mapping31 = Double.NaN;
			if(mappingPerc29.containsKey(sample))
			{
				mapping29 = mappingPerc29.get(sample)[0];
				outLine += mapping29+"\t";	 
			}
			else
				outLine += "-\t";
			if(mappingPerc31.containsKey(sample))
			{
				mapping31 = mappingPerc31.get(sample)[0];
				outLine += mapping31+"\t";
			}
			else
				outLine += "-\t";
			if(!Double.isNaN(mapping29) && !Double.isNaN(mapping31))
			{
				outLine += mapping29/mapping31+"\t";
			}
			else
				outLine += "-\t";
			if(mappingPerc31.containsKey(sample))
				outLine += mappingPerc31.get(sample)[1]+"\t";
			else if(mappingPerc29.containsKey(sample))
				outLine += mappingPerc29.get(sample)[1]+"\t";
			else
				outLine += "-\t";
			if(mappingPerc29.containsKey(sample))
				outLine += mappingPerc29.get(sample)[2]+"\t";
			else
				outLine += "-\t";
			if(mappingPerc31.containsKey(sample))
				outLine += mappingPerc31.get(sample)[2]+"\t";
			else
				outLine += "-\t";
			if(missing.containsKey(sample))
				outLine += missing.get(sample)+"\t";
			else
				outLine += "\t";
			if(failedKallisto.containsKey(sample))
				outLine += failedKallisto.get(sample)+"\t";
			else
				outLine += "\t";
			outLine+=line;
			writer.write(outLine+"\n");
		}
		reader.close();
		writer.close();
		System.out.println("Done! File written to: " + writeFN);
	}

	private static void elseif(boolean containsKey) {
		// TODO Auto-generated method stub
		
	}
}
