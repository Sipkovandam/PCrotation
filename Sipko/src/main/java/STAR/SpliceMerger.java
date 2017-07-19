package STAR;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.math3.stat.descriptive.rank.Median;

import Tools.FileSearcher;
import Tools.FileUtils;
import Tools.Script;

public class SpliceMerger extends Script<SpliceMerger>
{
	private String spliceFilesFolderComment = "/root/folder/results/; MANDATORY // Folder containing all the splice files from the 2nd pass";
	private String splicesFilesFolder = null;
	private String spliceFileExtencionComment = "SJ.out.tab; MANDATORY // string to search for when searching for splice files in the <resultsFolder>";
	private String spliceFileExtencion = "SJ.out.tab";
	private String writeFnComment = "/root/folder/SJ_Merged_2ndPass_annotated.out.tab; OPTIONAL// filename of the file in which the expression for all samples for all annotated splice junctions from the second pass should be written. First 6 columns of the splice file are collapsed into the first column becoming the names of the splice junctinos";
	private String writeFn = null;
	private String writeFN_Splice_2ndPassInputComment = "/root/directory/mergedSpliceSites.sj.out; optional // filenames of the input file for the splice junctions to be used in the second pass. These are all junctions in the fused in the first pass plus those where reads overlapping <= <readCutoff> reads";
	private String writeFN_Splice_2ndPassInput = null;// file where all the splice variants merged into 1 file with info on numbe of spliced reads overlapping each in all samples together
	
	private String writeSpliceSummaryFnComment = "/root/folder/SJ_Merged_2ndPass_annotated.out.tab; OPTIONAL// filename where the sum of read counts for each splice junction is written. This is written in a format that can be used by STAR for the 2nd pass";
	private String writeSpliceSummaryFn = null;
	private String annotatedOnlyComment = "true; if true includes only annotated junctions in the merged spliceFile";
	private boolean annotatedOnly = true;
	private String exonsInsteadOfJunctionsComment = "false; if true adapts to the featurecounts output file format (basically skips first 2 lines)";
	private boolean exonsInsteadOfJunctions = false;

	private String readCutoffComment = "8; The minimum number of reads a splice junction must be observed in, in order to be included";
	private int readSumCutoff = 8;
	
	private String sampleComment = "0; The minimum number of samples a splice junction must be observed in (with at least 1 read), in order to be included";
	private int sampleCutoff = 0;
	
	public void run()
	{
		try
		{
			FileSearcher fileSearcher = new FileSearcher();
			fileSearcher.setFolders(splicesFilesFolder);
			fileSearcher.setSearchStrings(new String[]{spliceFileExtencion});
			
			if(exonsInsteadOfJunctions)
				fileSearcher.setForbiddenStrings(new String[]{".out.summary"});
			
			if(writeFn == null)
				writeFn = new File(new File(splicesFilesFolder).getParent())+"/SJ_Merged.out.tab";
			
			if(writeFN_Splice_2ndPassInput == null)
				writeFN_Splice_2ndPassInput = new File(new File(splicesFilesFolder).getParent())+"/SJ_Merged_AnnotatedOrAboveReadCutoff.out.tab";
			
			if(writeSpliceSummaryFn == null)
				writeSpliceSummaryFn = new File(new File(splicesFilesFolder).getParent())+"/SJ_Merged_AnnotatedSummary.out.tab";
			
			String spliceFilesFn=new File(writeFn).getParent()+"/spliceFiles.txt";
			fileSearcher.setWriteName(spliceFilesFn);
			fileSearcher.run();
			
			int nSpliceFiles = countSpliceFiles(spliceFilesFn);
			
			//create a hash that contains the splice counts for aech sample for each junction
			BufferedReader spliceFnsReader = FileUtils.createReader(spliceFilesFn);
			
			String spliceFn = null;
			HashMap<String, int[]> spliceSiteToCounts = new HashMap<String, int[]>();
			
			int fileIndex = 0;// index of this file, will become the index of column in final mergedsplice-output file
			String[] spliceHeaders = new String[nSpliceFiles];
			
			HashMap<String, Boolean> isAnnotated = new HashMap<String, Boolean>(); 
			while((spliceFn=spliceFnsReader.readLine())!= null)
			{
				//get the splice file headern (folderName where the file is located)
				spliceHeaders[fileIndex]=new File(new File(spliceFn).getParent()).getName();
				//get the splice counts per splice junction (annotated junctions only)
				List<Object> returnObjects = addSpliceCounts(spliceFn, spliceSiteToCounts, nSpliceFiles, fileIndex, isAnnotated);
				spliceSiteToCounts = (HashMap<String, int[]>) returnObjects.get(0);
				isAnnotated = (HashMap<String, Boolean>) returnObjects.get(1);
				fileIndex++;
			}			
			
			writeSpliceSiteCountsPersSample(nSpliceFiles, spliceHeaders, spliceSiteToCounts, writeFn);
			writeSpliceSiteCountsSummaryAnnotatedCutoff(nSpliceFiles, spliceHeaders, spliceSiteToCounts, writeFN_Splice_2ndPassInput, isAnnotated);
			writeSpliceSiteCountsSummary(nSpliceFiles, spliceHeaders, spliceSiteToCounts, writeSpliceSummaryFn);
			
		}catch(Exception e){e.printStackTrace();}
		p("Done! File written to:" + writeFn);
	}
	
	private void writeSpliceSiteCountsSummaryAnnotatedCutoff(int nSpliceFiles,
																String[] spliceHeaders,
																HashMap<String, int[]> spliceSiteToCounts, 
																String writeFN_Splice_2ndPassInput2, 
																HashMap<String, Boolean> isAnnotated) throws FileNotFoundException, IOException
	{
		BufferedWriter spliceWriter = FileUtils.createWriter(writeFN_Splice_2ndPassInput);
		
		if(exonsInsteadOfJunctions)
		{
			String spliceHeader="SpliceSite\tSummedCounts\tExpressedInNSamples\tExonLength\tExpressed/ExonLength";
			spliceWriter.write(spliceHeader);
			spliceWriter.write("\n");
		}
		
		//write each splice site and corresponding values
		spliceSiteToCounts.forEach((spliceSite,
									values) ->
		{
			try
			{
				int sum = 0;
				int expressedSamples = 0;
				int nZeros=0;
				for (int val : values)
				{
					sum+=val;
					if(val>0)
						expressedSamples++;
					if(val==0)
						nZeros++;
				}
				
				if(sum>=readSumCutoff && expressedSamples>=sampleCutoff | (!exonsInsteadOfJunctions && isAnnotated.containsKey(spliceSite) && isAnnotated.get(spliceSite).equals("1")))
				{
					//write the splice name in the STAR format (to different columns)
					String[] spliceFeatures = spliceSite.split("_");
					if(exonsInsteadOfJunctions)
						spliceWriter.write(spliceSite+"\t");
					else
					{
						for(String spliceFeature:spliceFeatures)
						{
							spliceWriter.write(spliceFeature+"\t");
						}
					}
					
					//write the values
					spliceWriter.write(sum+"\t"+expressedSamples+"\t"+nZeros);

					if(exonsInsteadOfJunctions)//write length of exon in extra column
						spliceWriter.write("\t"+spliceFeatures[6]+"\t"+((double)sum)/Double.parseDouble(spliceFeatures[6]));
					spliceWriter.write("\n");
				}							
			} catch (Exception e)
			{
				e.printStackTrace();
			}
		});
		spliceWriter.close();
	}

	private void writeSpliceSiteCountsSummary(int nSpliceFiles, String[] spliceHeaders, HashMap<String, int[]> spliceSiteToCounts, String writeSpliceSummaryFn2) throws FileNotFoundException, IOException
	{
		p(""+writeSpliceSummaryFn);
		BufferedWriter spliceWriter = FileUtils.createWriter(writeSpliceSummaryFn);
		//write headers
		if(exonsInsteadOfJunctions)
		{
			String spliceHeader="SpliceSite\tSummedCounts\tExpressedInNSamples\tExonLength\tExpressed/ExonLength";
			spliceWriter.write(spliceHeader);
			spliceWriter.write("\n");
		}

		//write each splice site and corresponding values
		spliceSiteToCounts.forEach((spliceSite,
									values) ->
		{
			try
			{
				//write the splice name in the STAR format (to different columns)
				String[] spliceFeatures = spliceSite.split("_");
				if(exonsInsteadOfJunctions)
					spliceWriter.write(spliceSite+"\t");
				else
				{
					for(String spliceFeature:spliceFeatures)
					{
						spliceWriter.write(spliceFeature+"\t");
					}
				}
				
				//write the values
				int sum = 0;
				int expressedSamples = 0;
				int nZeros=0;
				for (int val : values)
				{
					sum+=val;
					if(val>0)
						expressedSamples++;
					if(val==0)
						nZeros++;
				}
				spliceWriter.write(sum+"\t"+expressedSamples+"\t"+ nZeros);
				if(exonsInsteadOfJunctions)//write length of exon in extra column
					spliceWriter.write("\t"+spliceFeatures[6]+"\t"+((double)sum)/Double.parseDouble(spliceFeatures[6]));
				spliceWriter.write("\n");
						
			} catch (Exception e)
			{
				e.printStackTrace();
			}
		});
		spliceWriter.close();
	}
	
	
	private void writeSpliceSiteCountsPersSample(int nSpliceFiles, String[] spliceHeaders, HashMap<String, int[]> spliceSiteToCounts, String writeFn2) throws FileNotFoundException, IOException
	{
		BufferedWriter spliceWriter = FileUtils.createWriter(writeFn);
		//write headers
		for (String spliceHeader : spliceHeaders)
		{
			spliceWriter.write("\t" + spliceHeader);
		}
		spliceWriter.write("\n");

		//write each splice site and corresponding values
		spliceSiteToCounts.forEach((k,
									v) ->
		{
			try
			{
				//write the splice name
				spliceWriter.write(k);
				//write the values
				for (int val : v)
				{
					spliceWriter.write("\t" + val);
				}
				spliceWriter.write("\n");
			} catch (Exception e)
			{
				e.printStackTrace();
			}
		});
		spliceWriter.close();
	}

	public String getspliceFilesFolder()
	{
		return splicesFilesFolder;
	}

	public void setspliceFilesFolder(String resultsFolder)
	{
		this.splicesFilesFolder = resultsFolder;
	}

	public String getSpliceFileExtencion()
	{
		return spliceFileExtencion;
	}

	public void setSpliceFileExtencion(String spliceFileExtencion)
	{
		this.spliceFileExtencion = spliceFileExtencion;
	}

	public String getWriteFn()
	{
		return writeFn;
	}

	public void setWriteFn(String writeFn)
	{
		this.writeFn = writeFn;
	}

	private List<Object> addSpliceCounts(	String spliceFn,
													HashMap<String, int[]> spliceSiteToCounts,
													int nSpliceFiles,
													int fileIndex, HashMap<String, Boolean> isAnnotated) throws IOException
	{
		BufferedReader spliceReader = FileUtils.createReader(spliceFn);
		String spliceLine = null;
		int annotationIndex = 5;
		int countIndex = 6;
		if(exonsInsteadOfJunctions)
		{
			spliceReader.readLine();spliceReader.readLine();
		}
		
		
		List<Object> returnObjects = null;
		while ((spliceLine = spliceReader.readLine()) != null)
		{
			returnObjects = addSpliceCount(spliceSiteToCounts,
												spliceLine,
												fileIndex,
												annotationIndex,
												countIndex,
												nSpliceFiles, isAnnotated);
		}

		return returnObjects;
	}

	private List<Object> addSpliceCount(	HashMap<String, int[]> spliceSiteToCounts,
													String spliceLine,
													int fileIndex,
													int annotationIndex,
													int countIndex,
													int nSpliceFiles, HashMap<String, Boolean> isAnnotated)
	{
		String[] features = spliceLine.split("\t");

		//if this feature does not have a non-ambiguous read overlapping it: skip it

		if (features[countIndex].equals("0"))
			return r(spliceSiteToCounts,isAnnotated);

		//if this feature is not annotated and these should be excluded: skip it
		if (!exonsInsteadOfJunctions && annotatedOnly && features[annotationIndex].equals("0"))
			return r(spliceSiteToCounts,isAnnotated);

		StringBuilder spliceNameBuilder = new StringBuilder();

		spliceNameBuilder.append(features[0]);
		if(exonsInsteadOfJunctions)//makes sure it is the same format as the splice junction files in which the gene names are also separated by a double _ from the rest
			spliceNameBuilder.append("_");
		for (int f = 1; f <= annotationIndex; f++)
		{
			spliceNameBuilder.append("_");
			spliceNameBuilder.append(features[f]);
		}
		String spliceName = spliceNameBuilder.toString();
		isAnnotated.put(spliceName, features[countIndex].equals("1"));
		
		int spliceCount = Integer.parseInt(features[countIndex]);

		//add the number to the array
		int[] countsArray = spliceSiteToCounts.get(spliceName);
		//if no array exists yet, create one
		if (countsArray == null)
			countsArray = new int[nSpliceFiles];

		countsArray[fileIndex] = spliceCount;
		//put it in the hash that contains all splice values
		spliceSiteToCounts.put(	spliceName,
								countsArray);

		return r(spliceSiteToCounts,isAnnotated);
	}

	private int countSpliceFiles(String spliceFilesFn) throws FileNotFoundException, IOException
	{
		int n = 0;
		BufferedReader spliceFnsReader = FileUtils.createReader(spliceFilesFn);
		while (spliceFnsReader.readLine() != null)
		{
			n++;
		}

		return n;
	}

	public String getSplicesFilesFolder()
	{
		return splicesFilesFolder;
	}

	public void setSplicesFilesFolder(String splicesFilesFolder)
	{
		this.splicesFilesFolder = splicesFilesFolder;
	}

	public boolean isAnnotatedOnly()
	{
		return annotatedOnly;
	}

	public void setAnnotatedOnly(boolean annotatedOnly)
	{
		this.annotatedOnly = annotatedOnly;
	}

	public String getWriteSpliceSummaryFn()
	{
		return writeSpliceSummaryFn;
	}

	public void setWriteSpliceSummaryFn(String writeSpliceSummaryFn)
	{
		this.writeSpliceSummaryFn = writeSpliceSummaryFn;
	}

	public int getReadSumCutoff()
	{
		return readSumCutoff;
	}

	public void setReadSumCutoff(int readSumCutoff)
	{
		this.readSumCutoff = readSumCutoff;
	}

	public int getSampleCutoff()
	{
		return sampleCutoff;
	}

	public void setSampleCutoff(int sampleCutoff)
	{
		this.sampleCutoff = sampleCutoff;
	}

	public String getWriteFN_Splice_2ndPassInput()
	{
		return writeFN_Splice_2ndPassInput;
	}

	public void setWriteFn_Splice_2ndPassInput(String writeFN_Splice_2ndPassInput)
	{
		this.writeFN_Splice_2ndPassInput = writeFN_Splice_2ndPassInput;
	}

	public boolean isExonsInsteadOfJunctions()
	{
		return exonsInsteadOfJunctions;
	}

	public void setExonsInsteadOfJunctions(boolean exonsInsteadOfJunctions)
	{
		this.exonsInsteadOfJunctions = exonsInsteadOfJunctions;
	}
}
