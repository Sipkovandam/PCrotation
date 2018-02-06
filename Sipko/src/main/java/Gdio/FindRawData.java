package Gdio;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Set;

import Slurm.ClusterHandler;
import Tools.FileSearcher;
import Tools.FileUtils;
import Tools.Script;

public class FindRawData extends Script<FindRawData>
{
	
	FileSearcher csvFileSearcher = new FileSearcher();
	FileSearcher rawDataFileSearcher = new FileSearcher();
	String writeFn = "test.txt";
	
	@Override
	public void run()
	{
		try
		{
			///Get all file locations
			log("Searching all raw data file locations");
			rawDataFileSearcher.run();//Need to uncomment this on the server
			
			log("Searching all samplesheet file locations");
			csvFileSearcher.run();
			
			log("Creating an overview of all sampleNames and corresponding file locations");
			//create hashtable <sampleName,paths(forward,backward)>
			HashMap<String, ArrayList<String>> fastqFile_To_FastqFileWithPath = readSampleNameToPathsHash(rawDataFileSearcher.getWriteName());
			
			//go through all samplesheet files
			String csvFilenamesFile = csvFileSearcher.getWriteName();
			BufferedReader csvFilenamesReader = FileUtils.createReader(csvFilenamesFile);
			String csvFn = null;

			BufferedWriter overviewWriter = FileUtils.createWriter(writeFn);
			overviewWriter.write("sampleName\tdnaNumber\tfastqLocationForward\tfastqLocationBackward\n");
			while((csvFn=csvFilenamesReader.readLine())!=null)
			{
				//for each samplesheet find the identifier and add it to the overview with the fastq files in the connected columns
				parseCsvAndCreateOverView(csvFn,fastqFile_To_FastqFileWithPath, overviewWriter);	
			}
			overviewWriter.close();
			log("Done, file written to: " + writeFn);
		}catch(Exception e){e.printStackTrace();}
			
	}

	//add all fastq files per samplename
	private HashMap<String, ArrayList<String>> readSampleNameToPathsHash(String writeName) throws IOException
	{
		BufferedReader fastqPathsReader = FileUtils.createReader(writeName);
		String fastqPath = null;
		HashMap<String,ArrayList<String>> fastqFile_To_FastqFilesWithPath= new HashMap<String,ArrayList<String>>();
		while((fastqPath=fastqPathsReader.readLine())!=null)
		{
			// /groups/umcg-gd/prm02/rawdata/ngs/140509_SN163_0551_BH8AN0ADXX/140509_SN163_0551_BH8AN0ADXX_L2_2.fq.gz

			//String sampleName = fastqPath.replaceAll("_.*_.*_*\\.fq\\.gz", "");
			String sampleName = fastqPath.replaceAll(".*/ngs/", "").replaceAll("_[1-2].fq.gz", ".fq.gz");
			log("sampleName=\t" + sampleName);
			ArrayList<String> fastqFiles = fastqFile_To_FastqFilesWithPath.get(sampleName);
			if(fastqFiles==null)
				fastqFiles=new ArrayList<String>();
			fastqFiles.add(fastqPath);
			fastqFile_To_FastqFilesWithPath.put(sampleName,fastqFiles);
		}
		return fastqFile_To_FastqFilesWithPath;
	}

	private void parseCsvAndCreateOverView(String csvFn, HashMap<String, ArrayList<String>> fastQfile_To_FastqFileWithPath, BufferedWriter overviewWriter)
	{
		try
		{
			BufferedReader csvReader = FileUtils.createReader(csvFn);
			
			//make a hash with the column names to column numbers
			String line = csvReader.readLine();
			HashMap<String, Integer> headerToColnumber =  FileUtils.makeHash(line, ",");
			
			HashMap<String,String> internalSampleId_To_ForwardFastqFiles = new HashMap<String,String>();
			HashMap<String,String> internalSampleId_To_BackwardFastqFiles = new HashMap<String,String>();
			HashMap<String, String> internalSampleId_To_DnaNumber = new HashMap<String, String>();
			while((line = csvReader.readLine())!=null)
			{
				//get all important features
				String[] eles = line.split(",");
				String internalSampleID = getColNumber("internalSampleID", headerToColnumber, eles, csvFn);
				String sequencingStartDate = getColNumber("sequencingStartDate", headerToColnumber, eles, csvFn);
				String sequencer = getColNumber("sequencer", headerToColnumber, eles, csvFn);
				String run = getColNumber("run", headerToColnumber, eles, csvFn);
				String flowcell = getColNumber("flowcell", headerToColnumber, eles, csvFn);
				String lane = "L"+getColNumber("lane", headerToColnumber, eles, csvFn);
				String barcode = getColNumber("barcode", headerToColnumber, eles, csvFn);
				String gender = getColNumber("Gender", headerToColnumber, eles, csvFn);
				String contact = getColNumber("contact", headerToColnumber, eles, csvFn);
				String project = getColNumber("project", headerToColnumber, eles, csvFn);
				
				if(internalSampleID==null)
				{
					log("internalSampleID missing for sample:\t" + line + "\n in file: \t" + csvFn);
					continue;
				}
				
				//get the DNA number 
				String dnaNumber = getDnaNumber(internalSampleID);
				
				internalSampleId_To_DnaNumber.put(internalSampleID, dnaNumber);

				//get the samplename which should be the foldername in the "/projects/rawData/ngs" folder
				//first the folder
				String runName = addToSampleName(sequencingStartDate, sequencer, run, flowcell);
				//then the fastqFileName
				String fastqName = addToSampleName(runName+ "/"+ runName, lane, barcode);
				fastqName+= ".fq.gz";
				
				//write overview
				ArrayList<String> fastqFilesWithPath = fastQfile_To_FastqFileWithPath.get(fastqName);
				
				if(!fastqFileExists(fastqFilesWithPath, internalSampleID, fastqName))
					continue;
				
				addToPerSampleHashes(internalSampleID, fastqFilesWithPath, internalSampleId_To_ForwardFastqFiles, internalSampleId_To_BackwardFastqFiles);

			}
			
			//create sample overview
			Set<String> internalSampleIds = internalSampleId_To_ForwardFastqFiles.keySet();
			for(String internalSampleId : internalSampleIds)
			{
				overviewWriter.write(internalSampleId+"\t" + internalSampleId_To_DnaNumber.get(internalSampleId) + "\t" +internalSampleId_To_ForwardFastqFiles.get(internalSampleId)+"\t"+internalSampleId_To_BackwardFastqFiles.get(internalSampleId)+"\n");
			}
			
			
			csvReader.close();
		} catch (FileNotFoundException e)
		{e.printStackTrace();} catch (IOException e)
		{e.printStackTrace();}
		
	}

	private void addToPerSampleHashes(	String internalSampleID, ArrayList<String> fastqFilesWithPath, HashMap<String, String> internalSampleId_To_ForwardFastqFiles,
										HashMap<String, String> internalSampleId_To_BackwardFastqFiles)
	{
		for(String fastqFileWithPath : fastqFilesWithPath)
			if(!fastqFileWithPath.endsWith("_2.fq.gz"))
				addToSampleOverviewHash(internalSampleID, fastqFileWithPath, internalSampleId_To_ForwardFastqFiles);
			else
			{
				addToSampleOverviewHash(internalSampleID, fastqFileWithPath, internalSampleId_To_BackwardFastqFiles);
			}
	}

	private void addToSampleOverviewHash(String internalSampleID, String fastqFileWithPath, HashMap<String, String> internalSampleId_To_FastqFiles)
	{
		String forwardFastqs = internalSampleId_To_FastqFiles.get(internalSampleID);
		if(forwardFastqs==null)
			forwardFastqs=fastqFileWithPath;
		else
			forwardFastqs+=","+fastqFileWithPath;
		internalSampleId_To_FastqFiles.put(internalSampleID, forwardFastqs);
	}

	private boolean fastqFileExists(ArrayList<String> fastqFile, String internalSampleID, String fastqName)
	{
		if(fastqFile == null)
		{
			log("No fastq files found for sample: " + internalSampleID + "\t" + fastqName);
			return false;
		}
		return true;
	}

	private String getColNumber(String colName, HashMap<String, Integer> headerToColnumber, String[] eles, String csvFn)
	{
		if(headerToColnumber.get(colName)==null)	
			return null;
		if(headerToColnumber.get(colName)>= eles.length)
		{
			log("Incorrect CSV file:\t" + csvFn);
			return null;
		}
		return eles[headerToColnumber.get(colName)];
	}

	private String addToSampleName(String... features)
	{
		//features are the bits that make up the sample name (like the sample date,run number, bar code, etc...)
		String sampleName = features[0];
		for(int f = 1; f < features.length; f++)
		{
			if(features[f]!=null && !features[f].equals("None"))
				sampleName+="_"+features[f];
		}
		return sampleName;
	}

	private String getDnaNumber(String internalSampleID)
	{
		// sed 's/.*_DNA/DNA/g' |  sed 's/_.*//g')
		if(internalSampleID==null)
			return null;
		String dnaNumber = internalSampleID.replaceAll(".*_DNA", "DNA").replaceAll("_.*", "");
		return dnaNumber;
	}
}



