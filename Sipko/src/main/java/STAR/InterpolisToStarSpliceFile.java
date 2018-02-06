package STAR;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;

import Tools.FileUtils;
import Tools.Script;

public class InterpolisToStarSpliceFile extends Script<InterpolisToStarSpliceFile>
{
	String fn = "E:/Groningen/Data/Annotation/GRCh37/intropolis.v1.hg19_100lines.tsv.gz";
	String writeFn = "E:/Groningen/Data/Annotation/GRCh37/intropolis.v1.hg19_100lines_STARformatted.tsv";
	int minNumberOfReads = 1;//minimum number of spanning reads that need to support the junction in order to be included
	int minNumberOfSamples = 1000;//1000 was used in the original paper analysis leading to  56,865 new junctions (18.6%) //minimum number of samples the splice site has to be detected in in order to be included
	
	//converts the file obtained from the paper, where they identified splice variants in 20k samples, so that it has the STAR format
	public void run() 
	{
		try
		{
			//String fn = "/groups/umcg-gdio/tmp04/umcg-svandam/Data/RNAseq/Annotation/GCRh37/intropolis.v1.hg19.tsv.gz";
			//String writeFn = "/groups/umcg-gdio/tmp04/umcg-svandam/Data/RNAseq/Annotation/GCRh37/intropolis.v1.hg19_STARformatted.tsv.gz";
			String writeFnExtra = FileUtils.removeExtention(writeFn) + "_extra.tsv.gz";		
			
			BufferedReader reader;
			BufferedWriter writer = FileUtils.createWriter(writeFn);
			BufferedWriter writerExtra = FileUtils.createWriter(writeFnExtra);
			reader = FileUtils.createReader(fn);
		
			String line = null;
			int n = 0;
			while((line=reader.readLine())!=null)
			{
				if(n<3)
					log("First 3 input lines :" + line);
				if(n%100000 == 0)
					log(n + " lines done");
				String[] eles = line.split("\t");
				String[] counts= eles[7].split(",");
				int spanningReads = 0;
				for(String count : counts)
					spanningReads += Integer.parseInt(count);
				
				String outline = eles[0].replace("chrM", "MT").replace("chr", "")+"\t"+eles[1]+"\t"+eles[2];
				if(eles[3].contains("+"))
						outline+="\t"+1;
				else if(eles[3].contains("-"))
						outline+="\t"+2;
				//GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5:AT/AC, 6: GT/AT
				if(eles[4].contains("GT") && eles[5].contains("AG"))
					outline+="\t"+1;
				if(eles[4].contains("CT") && eles[5].contains("AC"))
					outline+="\t"+2;
				if(eles[4].contains("GC") && eles[5].contains("AG"))
					outline+="\t"+3;
				if(eles[4].contains("CT") && eles[5].contains("GC"))
					outline+="\t"+4;
				if(eles[4].contains("AT") && eles[5].contains("AC"))
					outline+="\t"+5;
				if(eles[4].contains("GT") && eles[5].contains("AT"))
					outline+="\t"+6;
				outline+="\t"+1;//it is annotated
				outline+="\t"+spanningReads;//number of uniquel mapping reads
				outline+="\t"+0;//number of multimapping reads
				outline+="\t"+25;//overhang just set it to 25 as the other file lacks this information.
				int nSamples = eles[6].split(",").length;
				
				if(spanningReads>=minNumberOfReads && nSamples>=minNumberOfSamples)
				{	
					writer.write(outline+"\n");
					writerExtra.write(outline+"\t"+ nSamples+"\n");
				}	
				n++;
			}
			writer.close();
			writerExtra.close();
			reader.close();	
		}catch(Exception e){e.printStackTrace();}
	}
}

