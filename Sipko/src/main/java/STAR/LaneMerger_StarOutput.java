package STAR;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashMap;

import com.google.common.io.Files;

import Tools.FileUtils;
import Tools.SampleSheetGafParser;
import Tools.Script;

public class LaneMerger_StarOutput extends Script<LaneMerger_StarOutput>
{
	String starResultsFolderName = null;
	double mappingPercentageCutoff = 1;//currently not used. Should add somewhere in addOltToNewFileNameHash().
	String unMergedFolderName = null;//folder that will end up containing the folders when they were unmerged
	private String sampleSheetFn = null;
	
	@Override
	public void run()
	{
		try
		{
			init();
			
			SampleSheetGafParser gafSheet = new SampleSheetGafParser();
			gafSheet.setSampleSheetFn(sampleSheetFn);
			gafSheet.setSeparator(",");
			gafSheet.run();
			
			File[] files = new File(starResultsFolderName).listFiles();
			HashMap<String,ArrayList<String>> newFilesToOldFiles = addOldToNewFileNamesHash(files, gafSheet);
			
			p("newFilesToOldFiles.size="+newFilesToOldFiles.size());
			combineMultiLaneSamples(newFilesToOldFiles);
		}catch(Exception e){e.printStackTrace();}
	}

	private void init()
	{
		if(unMergedFolderName == null)
			unMergedFolderName = FileUtils.makeFolderNameEndWithSlash(new File(starResultsFolderName).getParent())+"Results_beforeMergeFolders/";	
	}

	private void combineMultiLaneSamples(HashMap<String, ArrayList<String>> newFolderToOldFolders) throws FileNotFoundException, IOException
	{
		for(String mergeFolder : newFolderToOldFolders.keySet())
		{
			if(newFolderToOldFolders.get(mergeFolder).size()<2)
				continue;
			
			FileUtils.makeDir(mergeFolder);
			
			ArrayList<String> oldFolders = newFolderToOldFolders.get(mergeFolder);
			
			String oldFolder = oldFolders.get(0);
			String[] filesToMerge = new File(oldFolder).list();
			for(String fileToMerge: filesToMerge)
			{
				merge(fileToMerge, oldFolders, mergeFolder);
			}

			p("Moving files to " + unMergedFolderName);
			new File(unMergedFolderName).mkdirs();
			
			//move the old unmerged folders to a different directory
			for(String oldFol: oldFolders)
			{
				p("Moving director\t=" + oldFol);
				String newDirectory = FileUtils.makeFolderNameEndWithSlash(unMergedFolderName)+new File(oldFol).getName();
				org.apache.commons.io.FileUtils.moveDirectory(new File(oldFol), new File(newDirectory));
			}
		}
	}

	private void merge(	String fileToMerge,
						ArrayList<String> oldFolders, String mergeFolder) throws IOException
	{
		String mergeFn= FileUtils.makeFolderNameEndWithSlash(mergeFolder)+fileToMerge;
		switch(fileToMerge.replace(".gz", ""))
		{
			case "SJ.out.tab":
				mergeSJFile(oldFolders, fileToMerge, mergeFn);
				break;
			case "ReadsPerGene.out.tab":
				mergeReadsPerGeneFile(oldFolders, fileToMerge, mergeFn);
				break;
			case "Log.progress.out":
				mergeLogProgressFile(oldFolders, fileToMerge, mergeFn);
				break;
			case "Log.out":
				mergeLogFile(oldFolders, fileToMerge, mergeFn);
				break;
			case "Log.final.out":
				mergeLogFinalStarFile(oldFolders, fileToMerge, mergeFn);
				break;
			case "featureCounts.out.summary":
				mergeFeatureCountsSummaryFile(oldFolders, fileToMerge, mergeFn);
				break;
			case "featureCounts.out":
				mergeFeatureCountsFile(oldFolders, fileToMerge, mergeFn);
				break;
		}
		
	}

	private void mergeFeatureCountsFile(ArrayList<String> oldFolders, String fileToMerge, String mergeFn) throws IOException
	{
		StringBuilder writeLineBuilder = new StringBuilder();
		BufferedWriter mergeWriter = FileUtils.createWriter(mergeFn);
		
		BufferedReader[] oldReaders = new BufferedReader[oldFolders.size()];
		for(int r = 0; r <oldFolders.size();r++)
		{
			oldReaders[r]=FileUtils.createReader(FileUtils.makeFolderNameEndWithSlash(oldFolders.get(r))+fileToMerge);
			String line = oldReaders[r].readLine();
			if(r==0)
				mergeWriter.write(line+"\n");
			line = oldReaders[r].readLine();
			if(r==0)
				mergeWriter.write(line+"_LanesMerged"+"\n");
		}
		
		
		String line = "";
		out: while(line!=null)
		{
			int value = 0;
			writeLineBuilder.setLength(0);
			for(int r = 0; r < oldReaders.length;r++)
			{
				line = oldReaders[r].readLine();
				if(line == null)
					break out;
				String[] eles = line.split("\t");
				
				try{
					value += Integer.parseInt(eles[6]);
				}catch(Exception e){p("line =" + line);e.printStackTrace();}
				
				if(r==0)
				{
					for(int x = 0; x < 6; x++)
					{
						writeLineBuilder.append(eles[x]);
						writeLineBuilder.append("\t");
					}
				}
			}
			writeLineBuilder.append(value);

			writeLineBuilder.append("\n");
			mergeWriter.write(writeLineBuilder.toString());
		}
		
		for(int r = 0; r <oldFolders.size();r++)
		{
			oldReaders[r].close();
		}
		mergeWriter.close();		
	}

	private void mergeFeatureCountsSummaryFile(	ArrayList<String> oldFolders,
												String fileToMerge,
												String mergeFn) throws IOException
	{
		StringBuilder writeLineBuilder = new StringBuilder();
		BufferedWriter mergeWriter = FileUtils.createWriter(mergeFn);
		
		BufferedReader[] oldReaders = new BufferedReader[oldFolders.size()];
		for(int r = 0; r <oldFolders.size();r++)
		{
			oldReaders[r]=FileUtils.createReader(FileUtils.makeFolderNameEndWithSlash(oldFolders.get(r))+fileToMerge);
			String line = oldReaders[r].readLine();
			if(r==0)
				mergeWriter.write(line+"_LanesMerged"+"\n");
		}
		
		
		String line = "";
		out: while(line!=null)
		{
			int value = 0;
			writeLineBuilder.setLength(0);
			
			for(int r = 0; r < oldReaders.length;r++)
			{
				line = oldReaders[r].readLine();
				if(line==null)
					break out;
				String[] eles = line.split("\t");

				value += Integer.parseInt(eles[1]);
				
				if(r==0)
				{
					writeLineBuilder.append(eles[0]);
					writeLineBuilder.append("\t");
				}
			}
			writeLineBuilder.append(value);

			writeLineBuilder.append("\n");
			mergeWriter.write(writeLineBuilder.toString());
		}

		for(int r = 0; r <oldFolders.size();r++)
		{
			oldReaders[r].close();
		}
		mergeWriter.close();	
	}

	private void mergeLogFinalStarFile(ArrayList<String> oldFolders, String fileToMerge, String mergeFn) throws FileNotFoundException, IOException
	{
		double mappingSpeed = 0;
		double inputReads = 0;
		double avgInputreadLength = 0;
		double mappedReads = 0;
		double mappedLengthSum = 0;
		double spliceTotal = 0;
		double spliceAnnotated = 0;
		double spliceGT_AG = 0;
		double spliceGC_AG = 0;
		double spliceAT_AC = 0;
		double splice_NonCanonica = 0;
		double mismatchPerBase = 0;
		double deletionPerBase = 0;
		double deletionAverageLength = 0;
		double insertionPerBase = 0;
		double insertionAverageLength = 0;
		
		double multiMappingReads = 0;
		double  tooManyLocusMappingReads = 0;
		
		double unmappedToManyMismatches = 0;
		double unmappedToShort = 0;
		double unmappedOther = 0;
		
		double chimericReads = 0;
		
		for(String oldFolder: oldFolders)
		{
			String oldFn = FileUtils.makeFolderNameEndWithSlash(oldFolder)+fileToMerge;
			BufferedReader oldFileReader = FileUtils.createReader(oldFn);
			String line = null;
			
			double sampleInputReads = 0;
			double sampleMappedReads = 0;
			double sampleAvgMappedLength = 0;
			while((line=oldFileReader.readLine())!=null)
			{
				String[] eles = line.split("\t");
				if(eles.length<2)
					continue;

				switch(eles[0])
				{
					case "       Mapping speed, Million of reads per hour |":
						mappingSpeed += Double.parseDouble(eles[1]);
						break;
					case "                          Number of input reads |":
						sampleInputReads = Integer.parseInt(eles[1]);
						inputReads += sampleInputReads;
						break;
					case "                      Average input read length |":
						avgInputreadLength += Double.parseDouble(eles[1]);
						break;
					case "                   Uniquely mapped reads number |":
						sampleMappedReads= Double.parseDouble(eles[1]);
						mappedReads += sampleMappedReads;
						break;
					case "                          Average mapped length |":
						sampleAvgMappedLength = Double.parseDouble(eles[1]);
						mappedLengthSum += sampleAvgMappedLength;
						break;
					case "                       Number of splices: Total |":
						spliceTotal += Double.parseDouble(eles[1]);
						break;
					case "            Number of splices: Annotated (sjdb) |":
						spliceAnnotated += Double.parseDouble(eles[1]);
						break;
					case "                       Number of splices: GT/AG |":
						spliceGT_AG += Double.parseDouble(eles[1]);
						break;
					case "                       Number of splices: GC/AG |":
						spliceGC_AG += Double.parseDouble(eles[1]);
						break;
					case "                       Number of splices: AT/AC |":
						spliceAT_AC += Double.parseDouble(eles[1]);
						break;
					case "               Number of splices: Non-canonical |":
						splice_NonCanonica += Double.parseDouble(eles[1]);
						break;
					case "                      Mismatch rate per base, % |":
						mismatchPerBase += Double.parseDouble(eles[1].replace("%", ""))*sampleMappedReads*sampleAvgMappedLength;
						break;
					case "                         Deletion rate per base |":
						deletionPerBase += Double.parseDouble(eles[1].replace("%", ""))*sampleMappedReads*sampleAvgMappedLength;
						break;
					case "                        Deletion average length |":
						deletionAverageLength += Double.parseDouble(eles[1])*sampleMappedReads;
						break;
					case "                        Insertion rate per base |":
						insertionPerBase += Double.parseDouble(eles[1].replace("%", ""))*sampleMappedReads*sampleAvgMappedLength;
						break;
					case "                       Insertion average length |":
						insertionAverageLength += Double.parseDouble(eles[1])*sampleMappedReads;
						break;
					case "        Number of reads mapped to multiple loci |":
						multiMappingReads += Double.parseDouble(eles[1]);
						break;
					case "        Number of reads mapped to too many loci |":
						tooManyLocusMappingReads += Double.parseDouble(eles[1]);
						break;
					case "       % of reads unmapped: too many mismatches |":
						unmappedToManyMismatches += Double.parseDouble(eles[1].replace("%", ""))*sampleInputReads;
						break;
					case "                 % of reads unmapped: too short |":
						unmappedToShort += Double.parseDouble(eles[1].replace("%", ""))*sampleInputReads;
						break;
					case "                     % of reads unmapped: other |":
						unmappedOther += Double.parseDouble(eles[1].replace("%", ""))*sampleInputReads;
						break;
					case "                       Number of chimeric reads |":
						chimericReads += Double.parseDouble(eles[1]);
						break;
				}
			}
			oldFileReader.close();
		}
		
		StringBuilder writeLineBuilder = new StringBuilder();
		BufferedWriter mergeWriter = FileUtils.createWriter(mergeFn);
		
		String oldFn = FileUtils.makeFolderNameEndWithSlash(oldFolders.get(0))+fileToMerge;
		BufferedReader oldFileReader = FileUtils.createReader(oldFn);
		String line = null;
		double nFiles = (double)oldFolders.size();
		while((line=oldFileReader.readLine())!=null)
		{
			writeLineBuilder.setLength(0);
			
			String[] eles = line.split("\t");
			writeLineBuilder.append(eles[0]);
			writeLineBuilder.append("\t");
			
			double mappedLength = mappedLengthSum/nFiles;
			
			switch(eles[0])
			{
				case "       Mapping speed, Million of reads per hour |":
					writeLineBuilder.append(mappingSpeed/nFiles);
					break;
				case "                          Number of input reads |":
					writeLineBuilder.append(inputReads);
					break;
				case "                      Average input read length |":
					writeLineBuilder.append(avgInputreadLength/nFiles);
					break;
				case "                   Uniquely mapped reads number |":
					writeLineBuilder.append(mappedReads);
					break;
				case "                        Uniquely mapped reads % |":
					writeLineBuilder.append(mappedReads/inputReads*100).append("%");
					break;
				case "                          Average mapped length |":
					writeLineBuilder.append(mappedLength);
					break;
				case "                       Number of splices: Total |":
					writeLineBuilder.append(spliceTotal);
					break;
				case "            Number of splices: Annotated (sjdb) |":
					writeLineBuilder.append((float)spliceAnnotated);
					break;
				case "                       Number of splices: GT/AG |":
					writeLineBuilder.append((float)spliceGT_AG);
					break;
				case "                       Number of splices: GC/AG |":
					writeLineBuilder.append((float)spliceGC_AG);
					break;
				case "                       Number of splices: AT/AC |":
					writeLineBuilder.append((float)spliceAT_AC);
					break;
				case "               Number of splices: Non-canonical |":
					writeLineBuilder.append((float) splice_NonCanonica );
					break;
				case "                      Mismatch rate per base, % |":
					writeLineBuilder.append(mismatchPerBase/mappedReads/mappedLength).append("%");
					break;
				case "                         Deletion rate per base |":
					writeLineBuilder.append(deletionPerBase/mappedReads/mappedLength).append("%");
					break;
				case "                        Deletion average length |":
					writeLineBuilder.append(deletionAverageLength/mappedReads);
					break;
				case "                        Insertion rate per base |":
					writeLineBuilder.append(insertionPerBase/mappedReads/mappedLength).append("%");
					break;
				case "                       Insertion average length |":
					writeLineBuilder.append(insertionAverageLength/mappedReads);
					break;
				case "        Number of reads mapped to multiple loci |":
					writeLineBuilder.append(multiMappingReads);
					break;
				case "             % of reads mapped to multiple loci |":
					writeLineBuilder.append(multiMappingReads/inputReads*100).append("%");
					break;
				case "        Number of reads mapped to too many loci |":
					writeLineBuilder.append(tooManyLocusMappingReads);
					break;
				case "             % of reads mapped to too many loci |":
					writeLineBuilder.append(tooManyLocusMappingReads/inputReads*100).append("%");
					break;
				case "       % of reads unmapped: too many mismatches |":
					writeLineBuilder.append(unmappedToManyMismatches/inputReads).append("%");
					break;
				case "                 % of reads unmapped: too short |":
					writeLineBuilder.append(unmappedToShort/inputReads).append("%");
					break;
				case "                     % of reads unmapped: other |":
					writeLineBuilder.append(unmappedOther/inputReads).append("%");
					break;
				case "                       Number of chimeric reads |":
					writeLineBuilder.append(chimericReads);
					break;
				default:
					writeLineBuilder.setLength(0);
					writeLineBuilder.append(line);
					break;
			}
			writeLineBuilder.append("\n");
			mergeWriter.write(writeLineBuilder.toString());
			
		}
		oldFileReader.close();
		mergeWriter.close();
	}

	private void mergeLogFile(ArrayList<String> oldFolders, String fileToMerge, String mergeFn) throws IOException
	{
		String oldFn = FileUtils.makeFolderNameEndWithSlash(oldFolders.get(0))+fileToMerge;
		Files.copy(new File(oldFn), new File(mergeFn));
	}

	private void mergeLogProgressFile(ArrayList<String> oldFolders, String fileToMerge, String mergeFn) throws FileNotFoundException, IOException
	{
		BufferedWriter mergeWriter = FileUtils.createWriter(mergeFn);
		for(int f = 0; f <oldFolders.size();f++)
		{
			String oldFn = FileUtils.makeFolderNameEndWithSlash(oldFolders.get(f))+fileToMerge;
			BufferedReader oldFileReader = FileUtils.createReader(oldFn);
			String line = null;
			if(f != 0)
			{
				oldFileReader.readLine();
				oldFileReader.readLine();
			}
			
			while((line=oldFileReader.readLine())!=null)
			{
				
				if(line.startsWith("ALL DONE!") && f != oldFolders.size()-1)
					continue;
					
				mergeWriter.write(line+"\n");
			}
			oldFileReader.close();
		}
		mergeWriter.close();
	}

	private void mergeReadsPerGeneFile(ArrayList<String> oldFolders, String fileToMerge, String mergeFn) throws IOException
	{
		StringBuilder writeLineBuilder = new StringBuilder();
		BufferedWriter mergeWriter = FileUtils.createWriter(mergeFn);
		
		BufferedReader[] oldReaders = new BufferedReader[oldFolders.size()];
		for(int r = 0; r <oldFolders.size();r++)
		{
			oldReaders[r]=FileUtils.createReader(FileUtils.makeFolderNameEndWithSlash(oldFolders.get(r))+fileToMerge);
		}
		
		String geneName = null;
		
		String line = "";
		out: while(line!=null)
		{
			int[] values = new int[3]; 
			writeLineBuilder.setLength(0);
			for(int r = 0; r < oldReaders.length;r++)
			{
				line = oldReaders[r].readLine();
				if(line == null)
					break out;
				String[] eles = line.split("\t");
				geneName=eles[0];
				
				for(int e=1; e< eles.length;e++)
				{
					values[e-1] += Integer.parseInt(eles[e]);
				}
			}
			writeLineBuilder.append(geneName);
			for(int value: values)
			{
				writeLineBuilder.append("\t");
				writeLineBuilder.append(value);
			}
			writeLineBuilder.append("\n");
			mergeWriter.write(writeLineBuilder.toString());
		}
		
		for(int r = 0; r <oldFolders.size();r++)
		{
			oldReaders[r].close();
		}
		mergeWriter.close();
	}

	private void mergeSJFile(	ArrayList<String> oldFolders,
								String fileToMerge, String mergeFn) throws FileNotFoundException, IOException
	{
		HashMap<String,int[]> spliceToValues = new HashMap<String,int[]>();
		
		for(String oldFolder: oldFolders)
		{
			String oldFn = FileUtils.makeFolderNameEndWithSlash(oldFolder)+fileToMerge;
			BufferedReader oldFileReader = FileUtils.createReader(oldFn);
			String line = null;
			StringBuilder spliceNameBuilder = new StringBuilder();
			String spliceName = null;
			while((line=oldFileReader.readLine())!=null)
			{
				String[] eles = line.split("\t");

				spliceNameBuilder = getSpliceName(eles, spliceNameBuilder);
				spliceName =spliceNameBuilder.toString();
				int[] values = spliceToValues.get(spliceName);
				if(values == null)
					values = new int[3];
				values[0]+=Integer.parseInt(eles[6]);//reads overlapping this feature only
				values[1]+=Integer.parseInt(eles[7]);//reads overlapping multiple features
				int maxOverhang = Integer.parseInt(eles[8]);
				if(maxOverhang > values[2])
					values[2]=maxOverhang;
				spliceToValues.put(spliceName.toString(), values);
			}
			oldFileReader.close();
		}
		
		StringBuilder writeLineBuilder = new StringBuilder();
		BufferedWriter mergeWriter = FileUtils.createWriter(mergeFn);
		for(String spliceName:spliceToValues.keySet())
		{
			int[] values = spliceToValues.get(spliceName);
			writeLineBuilder.setLength(0);
			
			writeLineBuilder.append(spliceName);
			for(int value : values)
			{
				writeLineBuilder.append("\t");
				writeLineBuilder.append(value);
			}
			writeLineBuilder.append("\n");
			mergeWriter.write(writeLineBuilder.toString());
			
		}
		
		mergeWriter.close();
	}

	private StringBuilder getSpliceName(String[] eles,
										StringBuilder spliceName)
	{
		spliceName.setLength(0);
		spliceName.append(eles[0]);
		for(int e = 1; e < 6;e++)
		{
			spliceName.append("\t");
			spliceName.append(eles[e]);
		}
		
		return spliceName;
	}

	private boolean checkMappingPercentage(String oldName) throws FileNotFoundException, IOException
	{
		String mappingPercentageStr = getMappingPercentage(oldName);
		if(mappingPercentageStr==null)
			return false;
		double mappingPercentage = Double.parseDouble(mappingPercentageStr);
		if(mappingPercentage>mappingPercentageCutoff)
			return true;
		return false;
	}

	private String getMappingPercentage(String oldName) throws FileNotFoundException, IOException
	{
		BufferedReader reader = FileUtils.createReader(oldName);
		String line = null;
		while((line=reader.readLine())!=null)
		{
			if(line.startsWith("                        Uniquely mapped reads %"))
				continue;
			String readsStr = line.split("\t")[1].replace("%","");
			return readsStr;
		}
		p("Corrupt mapping file detected: " + oldName);
		return null;
	}

	private HashMap<String, ArrayList<String>> addOldToNewFileNamesHash(File[] files, SampleSheetGafParser gafSheet)
	{
		HashMap<String,ArrayList<String>> newFilesToOldFiles = new HashMap<String,ArrayList<String>>();
		HashMap<String,String> sumFolderNameToSampelId = new HashMap<String,String>();//this becomes a problem when there are multiple different samples that are run on multiple lanes but have the same barcode (Patrick's idea to drop the laneNumbers biting us in the ass here)
		for(File file:files)
		{
			String folderName = file.getAbsolutePath();//foldername == sampleId
			String sumFolderName = folderName.replaceAll("_L[0-9]{1,2}_", "_");//doesnt work if lane number become bigger then 100
			
			if(!isSameInternalId(sumFolderNameToSampelId, gafSheet, folderName, sumFolderName))
				continue;
			
			ArrayList<String> oldnames = newFilesToOldFiles.get(sumFolderName);
			oldnames=FileUtils.addStringToArrayList(oldnames, folderName);
			newFilesToOldFiles.put(sumFolderName,oldnames);
		}
		
		//remove the sumNames from any list that is larger than 1 (can happen when you run code 2X)
		for(String newName : newFilesToOldFiles.keySet())
		{
			ArrayList<String> oldFiles =  newFilesToOldFiles.get(newName);
			if(oldFiles.size()>1)
				oldFiles.remove(newName);
			newFilesToOldFiles.put(newName, oldFiles);
		}	
		
		return newFilesToOldFiles;
	}

	private boolean isSameInternalId(HashMap<String, String> sumFolderNameToSampelId, SampleSheetGafParser gafSheet, String folderName, String sumFolderName)
	{
		String internalId = gafSheet.getCol("internalSampleID").get(new File(folderName).getName());//get the internalId of this sample
		if(internalId ==null)
			return false;

		String internalIdOtherSample = sumFolderNameToSampelId.get(sumFolderName);
		if(internalIdOtherSample==null)
		{
			internalIdOtherSample=internalId;
			sumFolderNameToSampelId.put(sumFolderName, internalIdOtherSample);
		}

		return internalIdOtherSample.equals(internalId);
	}

	public String getStarResultsFolderName()
	{
		return starResultsFolderName;
	}

	public void setStarResultsFolderName(String starResultsFolderName)
	{
		this.starResultsFolderName = starResultsFolderName;
	}

	public String getUnMergedFolderName()
	{
		return unMergedFolderName;
	}

	public void setUnMergedFolderName(String unMergedFolderName)
	{
		this.unMergedFolderName = unMergedFolderName;
	}
}
