package TextEditing;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import MatrixScripts.MatrixString;
import Tools.FileSearcher;
import Tools.FileUtils;
import Tools.Script;

public class BiosSampleOverviewMaker extends Script<BiosSampleOverviewMaker>
{
	FileSearcher fileSearcher = new FileSearcher();
	String sampleLijst = "";
	String conversionFile = "";//some weird file with a column contain
	String writeFn= "";//where the conversion table is written (from runId to file location);
	String writeFnMissingSamples= "";//where the conversion table is written (from runId to file location);
	String forwardName="_R1.fq.gz";
	String backwardName="_R2.fq.gz";
	String ourSamples="/groups/umcg-gdio/tmp04/umcg-svandam/Data/RNAseq/FastqFiles/";
	
	@Override
	public void run()
	{
		try
		{
			//search all the files
			fileSearcher.run();
			
			HashMap<String, Integer> samplesToInclude = FileUtils.StringToIndexHash(sampleLijst,false);
			
			//create conversion hash frokm runId to lldeep number. (example runID: AD1H7HACXX-1-6). lldeep number example: 103001208027
			HashMap<String, ArrayList<String>> lldeepNumberToRunIds= FileUtils.readStringMultiStringArrayList(conversionFile,10,1,false);
			
			BufferedReader fnReader = FileUtils.createReader(fileSearcher.getWriteName());
			BufferedWriter conversionWriter = FileUtils.createWriter(writeFn);
			
			String fn_full = null;
			while((fn_full=fnReader.readLine())!=null)
			{
				String fn = new File(fn_full).getName();
				String[] fn_split = fn.split("_");
				
				if(fn_split.length<5 || fn_full.contains(ourSamples))
				{
					String sampleName = fn.replace(this.forwardName, "").replace(this.backwardName, "");
					if(samplesToInclude.containsKey(sampleName))
					{
						samplesToInclude.remove(sampleName);
						
						addPairedOrSingleEndSampleToFile(sampleName,fn_full, conversionWriter, this.forwardName, this.backwardName);
					}
				}
				else//if it has a lifelines deep number, do an extra conversion
				{
					//llNumberToSampleName.put(fn_split[4], fn);
					String llDeepNumber = fn_split[4];
					//p("llDeepNumber="  + llDeepNumber);
					ArrayList<String> runIds= lldeepNumberToRunIds.get(llDeepNumber);
					if(runIds==null)//maybe I should put a warning here?
					{
						//System.out.println("llDeepNumber: "+ llDeepNumber + " not found in conversion file");
						continue;
					}
					for(String runId:runIds)
					{
						if(samplesToInclude.containsKey(runId))
						{
							samplesToInclude.remove(runId);
							addPairedOrSingleEndSampleToFile(runId,fn_full, conversionWriter, this.forwardName, this.backwardName);
						}
					}	
				}
			}
			
			fnReader.close();
			conversionWriter.close();
			
			BufferedWriter missingSampleWriter = FileUtils.createWriter(this.writeFnMissingSamples); 
			for(String sample:samplesToInclude.keySet())
			{
				missingSampleWriter.write(sample+"\n");
			}
			missingSampleWriter.close();
			log("Done, file written to:\n" + writeFn);
		}catch(Exception e){e.printStackTrace();}
	}

	private void addPairedOrSingleEndSampleToFile(String sampleName, String fn_full, BufferedWriter conversionWriter, String forwardName, String backwardName) throws IOException
	{
		if(fn_full.contains(forwardName) || fn_full.contains(backwardName))
		{
			conversionWriter.write(sampleName+"\t"+fn_full.replace(backwardName, forwardName)+"\n");
			conversionWriter.write(sampleName+"\t"+fn_full.replace(forwardName, backwardName)+"\n");
		}
		else
			conversionWriter.write(sampleName+"\t"+fn_full+"\n");
	}
}
