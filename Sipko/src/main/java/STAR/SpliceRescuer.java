package STAR;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import Tools.FileSearcher;
import Tools.FileUtils;
import Tools.Script;

public class SpliceRescuer extends Script<SpliceRescuer>
{
	private static final long serialVersionUID = -5902739316302729000L;

	String studyFnsFn = null;//all files from which you would like to retrieve new splice sites
	FileSearcher spliceFileSearcher = new FileSearcher();//SJ.out.tab.gz (all files you want to inlcude in the output0
	int minCountForInclude = 8;//splice site needs to have at least this many reads in any of the study samples
	int maxSamplesExpressed = 1;
	String writeFn = null;

	@Override
	public void run()
	{
		try
		{
			spliceFileSearcher.run();

			if (writeFn == null)
				writeFn = FileUtils.removeExtention(studyFnsFn) + "_rescueSites.txt.gz";

			HashMap<String, HashMap<String, Integer>> spliceSitesToSamplesToCounts = new HashMap<String, HashMap<String, Integer>>();

			boolean addNewSpliceSites = true;
			addCountsToHash(studyFnsFn,
							spliceSitesToSamplesToCounts,
							addNewSpliceSites,
							minCountForInclude);

			addNewSpliceSites = false;

			addCountsToHash(spliceFileSearcher.getWriteName(),
							spliceSitesToSamplesToCounts,
							addNewSpliceSites);
			System.out.println("n spliceSitesToSamplesToCounts = " + spliceSitesToSamplesToCounts.size());
			writeSpliceSitesToRescue(	spliceSitesToSamplesToCounts,
										spliceFileSearcher.getWriteName());
			
			p("Output written at:\t" + writeFn);

		} catch (Exception e)
		{
			e.printStackTrace();
		}
	}

	private void writeSpliceSitesToRescue(	HashMap<String, HashMap<String, Integer>> spliceSitesToSamplesToCounts,
											String includedFns) throws IOException
	{
		HashMap<String, Integer> sampleNameToColumn = new HashMap<String, Integer>();
		BufferedWriter spliceWriter = writeHeaderAndAndSetSampleNamesToColumnHash(	includedFns,
																					sampleNameToColumn,
																					this.writeFn);

		for (String spliceSite : spliceSitesToSamplesToCounts.keySet())
		{
			StringBuilder spliceLine = new StringBuilder();
			HashMap<String, Integer> samplesToCounts = spliceSitesToSamplesToCounts.get(spliceSite);

			if (samplesToCounts.size() > maxSamplesExpressed)
				continue;

			spliceLine.append(spliceSite);

			spliceLine = addCounts(	spliceLine,
									sampleNameToColumn,
									samplesToCounts);

			spliceLine.append("\n");
			spliceWriter.write(spliceLine.toString());
		}
		spliceWriter.close();
	}

	private StringBuilder addCounts(StringBuilder spliceLine,
									HashMap<String, Integer> sampleNameToColumn,
									HashMap<String, Integer> samplesToCounts)
	{
		int[] rowCounts = new int[sampleNameToColumn.size()];
		for (String sampleName : sampleNameToColumn.keySet())
		{
			int col = sampleNameToColumn.get(sampleName);
			int count = 0;
			if (samplesToCounts.get(sampleName) != null)
				count = samplesToCounts.get(sampleName);
			rowCounts[col] = count;
		}

		for (int rowCount : rowCounts)
		{
			spliceLine.append("\t");
			spliceLine.append(rowCount);
		}
		return spliceLine;
	}

	private BufferedWriter writeHeaderAndAndSetSampleNamesToColumnHash(	String includedFns,
																		HashMap<String, Integer> sampleNameToColumn,
																		String writeFn) throws FileNotFoundException, IOException
	{
		BufferedWriter spliceWriter = FileUtils.createWriter(writeFn);
		BufferedReader includeFnsReader = FileUtils.createReader(includedFns);
		String fn = null;
		int c = 0;
		StringBuilder header = new StringBuilder();
		header.append("spliceName");
		while ((fn = includeFnsReader.readLine()) != null)
		{
			String sampleName = FileUtils.getFolderName(fn);
			sampleNameToColumn.put(	sampleName,
									c);
			header.append("\t");
			header.append(sampleName);
			c++;
		}
		header.append("\n");
		spliceWriter.write(header.toString());
		return spliceWriter;
	}

	private void addCountsToHash(	String studyFnsFn2,
									HashMap<String, HashMap<String, Integer>> spliceSitesToSamplesToCounts,
									boolean addNewSpliceSites) throws FileNotFoundException, IOException
	{
		addCountsToHash(studyFnsFn2,
						spliceSitesToSamplesToCounts,
						addNewSpliceSites,
						0);
	}

	private void addCountsToHash(	String studyFnsFn2,
									HashMap<String, HashMap<String, Integer>> spliceSitesToSamplesToCounts,
									boolean addNewSpliceSites,
									int countCutoff) throws FileNotFoundException, IOException
	{
		BufferedReader studyFnsReader = FileUtils.createReader(studyFnsFn);
		String studyFn = null;

		HashMap<String, ArrayList<String>> newFilesToOldFiles = addOldToNewFileNamesHash(studyFnsFn);

		while ((studyFn = studyFnsReader.readLine()) != null)
		{
			addSpliceSitesToHash(	studyFn,
									spliceSitesToSamplesToCounts,
									addNewSpliceSites,
									countCutoff,
									newFilesToOldFiles);
		}
	}

	private HashMap<String, ArrayList<String>> addOldToNewFileNamesHash(String studyFnsFn) throws IOException
	{
		HashMap<String, ArrayList<String>> newFilesToOldFiles = new HashMap<String, ArrayList<String>>();
		BufferedReader studyReader = FileUtils.createReader(studyFnsFn);
		String line = null;
		while ((line = studyReader.readLine()) != null)
		{
			String folderName = FileUtils.getFolderName(line);

			String sumFolderName = folderName.replaceAll("_L[0-9]{1,2}_", "_");//doesnt work if lane number become bigger then 100
			ArrayList<String> oldnames = newFilesToOldFiles.get(sumFolderName);
			FileUtils.addStringToArrayList(	oldnames,
											folderName);
		}
		return newFilesToOldFiles;
	}

	private void addSpliceSitesToHash(	String studyFn,
										HashMap<String, HashMap<String, Integer>> spliceSitesToSamplesToCounts,
										boolean addNewSites,
										int countCutoff,
										HashMap<String, ArrayList<String>> newFilesToOldFiles) throws FileNotFoundException, IOException
	{
		String spliceLine = null;
		BufferedReader[] spliceReaders = createReaders(studyFn, newFilesToOldFiles);
	
		String sampleName = FileUtils.getFolderName(studyFn);
		if(spliceReaders.length>1)
			sampleName=sampleName.replaceAll("_L[0-9]{1,2}_", "_");
		
		while ((spliceLine = spliceReaders[0].readLine()) != null)
		{
			String[] eles = spliceLine.split("\t");
			String annotated = eles[5];
			if (annotated.equals("1"))
				continue;

			int readCount = Integer.parseInt(eles[6]);
			readCount = addSpliceIntCountsFromOtherFiles(spliceReaders, readCount);
			
			if (readCount == 0)
				continue;

			String spliceSite = toSpliceSite(eles);
			HashMap<String, Integer> samplesToCounts = spliceSitesToSamplesToCounts.get(spliceSite);
			if (samplesToCounts == null && addNewSites == false)
				continue;

			if (samplesToCounts == null && readCount < countCutoff)
				continue;

			if (samplesToCounts == null)
				samplesToCounts = new HashMap<String, Integer>();

			samplesToCounts.put(sampleName,
								readCount);

			spliceSitesToSamplesToCounts.put(	spliceSite,
												samplesToCounts);
		}
	}

	private int addSpliceIntCountsFromOtherFiles(	BufferedReader[] spliceReaders,
													int readCount) throws IOException
	{
		if(spliceReaders.length>1)
		{
			for(int r =1; r< spliceReaders.length; r++)
			{
				String spliceLine = spliceReaders[r].readLine();
				String[] eles = spliceLine.split("\t");
				readCount += Integer.parseInt(eles[6]);
				System.out.println(readCount);
			}
		}

		return readCount;
	}

	private BufferedReader[] createReaders(String studyFn, HashMap<String, ArrayList<String>> newFilesToOldFiles) throws FileNotFoundException, IOException
	{

		String potentialNewName = studyFn.replaceAll("_L[0-9]{1,2}_", "_");
		ArrayList<String> oldFileNames = newFilesToOldFiles.get(potentialNewName);
		if(oldFileNames==null)
			return new BufferedReader[]{FileUtils.createReader(studyFn)};
		
		BufferedReader[] spliceReaders = new BufferedReader[oldFileNames.size()];
		for(int f = 0; f < oldFileNames.size(); f++)
			spliceReaders[f] = FileUtils.createReader(oldFileNames.get(f));
		return spliceReaders;
	}

	private String toSpliceSite(String[] eles)
	{
		StringBuilder spliceBuilder = new StringBuilder();
		spliceBuilder.append(eles[0]);
		for (int e = 1; e < 6; e++)
		{
			spliceBuilder.append("_");
			spliceBuilder.append(eles[e]);
		}
		return spliceBuilder.toString();
	}
}
