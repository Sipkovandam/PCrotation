package R;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import RowAnalyses.RowAboveCutoffCounter;
import RowAnalyses.RowAverageCalculator;
import RowAnalyses.RowBelowCutoffCounter;
import RowAnalyses.RowColumnGetter;
import RowAnalyses.RowJob;
import RowAnalyses.RowJobExecutor;
import RowAnalyses.RowMaxGetter;
import RowAnalyses.RowMedianGetter;
import RowAnalyses.RowStdevCalculator;
import RowAnalyses.RowZscoreCalculator;
import Tools.FileUtils;
import Tools.Script;
import Tools.StringFilter;
import Tools.StringFilters;
import Tools.StringFilters.ForbiddenStringChecker;
import Tools.StringFilters.RequiredStringChecker;

public class SplicingAnalysisPipeline2 extends Script<SplicingAnalysisPipeline>
{
	//creates files per patient with z-scores based on BIOS averages and stdevs and z-scores based on non-family members of the same study averages and stdevs.
	//Also adds gene network scores and mutations
	String correctedMatrixFn = null;
	String uncorrectedMatrixFn = null;
	String rawMatrixFn = null;

	String[] ourStudyStrings = new String[] { "SN163", "NB501043_0096" };

	String gnFolder = null;
	String patientSampleToFamilyIdFn = null;
	String parentSampleToFamilyIdFn = null;
	String ensgNameToGeneSymbolFn = null;
	String enstToEnsgFn = null;
	String vcfFolder = null;
	String sampleNameToFamilySampleNamesFn = null;//has patient sampleId in first column and sampleId (including patient) of all family members in second column (so 3 lines per patient)
	int summaryCutoff = 2; //any genes with a z-scores larger then this value will be in the summary file
	String writeFolder = null;
	double maxOutlierCount = 4; //maximum number of times a gene is allowed to be an outleir in the other parents. Also maximum number of times a gene is allowed to be an outleir in the other cases
	double medianCutoff = 8;
	int nThreads = 1;//Do not run this multithreaded. I got halfway implementing it, but did not finish. Some writeVariables are used by multiple threads simultaneously which is dangerous.
	private transient int zScoreCutoff = 3;//cutoff die Patrick gebruikt voor het excluderen van outliers die vaker voorkomen.

	@Override
	public void run()
	{
		try
		{
			HashMap<String, String> patientSampleNameToFamilyId = FileUtils.readStringStringHash(	patientSampleToFamilyIdFn,
																									1);
			HashMap<String, String> parentSampleNameToFamilyId = FileUtils.readStringStringHash(parentSampleToFamilyIdFn,
																								1);
			HashMap<String, HashMap<String, Integer>> caseSampleNameToFamilySampleNames = FileUtils.readStringMultiStringHash(	sampleNameToFamilySampleNamesFn,
																																0,
																																1,
																																false,
																																true);

			HashMap<String, String> ensgNameToGeneSymbol = FileUtils.readStringStringHash(	ensgNameToGeneSymbolFn,
																							false,
																							1,
																							1,
																							0);

			//get all the different index lists over which the statistics should be calculated
			StringFilters stringFiltersBios = new StringFilters();
			stringFiltersBios.addFilter(stringFiltersBios.new ForbiddenStringChecker(ourStudyStrings));
			int[] biosIndexes = getIndexes(	this.rawMatrixFn,
											stringFiltersBios);
			StringFilters stringFiltersStudy = new StringFilters();
			stringFiltersStudy.addFilter(stringFiltersStudy.new RequiredStringChecker(ourStudyStrings));
			int[] studyIndexes = getIndexes(this.rawMatrixFn,
											stringFiltersStudy);

			HashMap<String, int[]> sampleToOtherIndexesStudy = getStudyIndexesPerSample(this.rawMatrixFn,
																						stringFiltersStudy,
																						caseSampleNameToFamilySampleNames);

			createWriteFolder();

			String biosBasedResultsFolder = this.writeFolder + "bios/";

			//jobs to be executed on each row
			String averagesFn = "averages.txt";
			RowAverageCalculator rowAverageCalculator = new RowAverageCalculator();
			rowAverageCalculator.setWriteFn(averagesFn);

			String stDevFn = "standardDeviations.txt";
			RowStdevCalculator rowStdevCalculator = new RowStdevCalculator();
			rowStdevCalculator.setWriteFn(stDevFn);

			RowZscoreCalculator rowBiosBasedZscoreCalculator = new RowZscoreCalculator();
			String biosBasedZscoreFn = "biosBasedZscores.txt";
			rowBiosBasedZscoreCalculator.setWriteFn(biosBasedZscoreFn);
			ArrayList<RowJobExecutor> rowJobExecutors = new ArrayList<RowJobExecutor>();

			//bios indexes
			FileUtils.makeDir(biosBasedResultsFolder);
			RowJobExecutor rowJobExecutorBiosIndexes = new RowJobExecutor();
			rowJobExecutorBiosIndexes.addJob(rowAverageCalculator);
			rowJobExecutorBiosIndexes.addJob(rowStdevCalculator);
			rowJobExecutorBiosIndexes.addJob(rowBiosBasedZscoreCalculator);
			rowJobExecutorBiosIndexes.setIncludeIndexes(biosIndexes);
			rowJobExecutorBiosIndexes.setWriteFolder(biosBasedResultsFolder);

			///////////////////////////////////////////////////////////////////////////////////////////////////////////
			//Study indexes
			//jobs to be executed on each row
			RowAverageCalculator rowAverageCalculatorStudy = new RowAverageCalculator();
			rowAverageCalculatorStudy.setWriteFn(averagesFn);

			RowStdevCalculator rowStdevCalculatorStudy = new RowStdevCalculator();
			rowStdevCalculatorStudy.setWriteFn(stDevFn);

			RowZscoreCalculator rowBiosBasedZscoreCalculatorStudy = new RowZscoreCalculator();
			String biosBasedZscoreFnStudy = "biosBasedZscores.txt";
			rowBiosBasedZscoreCalculatorStudy.setWriteFn(biosBasedZscoreFnStudy);
			rowBiosBasedZscoreCalculatorStudy.setExecutorWithAvgStdevsToUse(rowJobExecutorBiosIndexes);
			rowBiosBasedZscoreCalculatorStudy.setRowAverageCalculator(rowAverageCalculatorStudy);
			rowBiosBasedZscoreCalculatorStudy.setRowStdevCalculator(rowStdevCalculatorStudy);

			String studyBasedResultsFolder = this.writeFolder + "study/";
			FileUtils.makeDir(studyBasedResultsFolder);
			RowJobExecutor rowJobExecutorStudyIndexes = new RowJobExecutor();
			rowJobExecutorStudyIndexes.addJob(rowAverageCalculatorStudy);
			rowJobExecutorStudyIndexes.addJob(rowStdevCalculatorStudy);
			rowJobExecutorStudyIndexes.addJob(rowBiosBasedZscoreCalculatorStudy);
			rowJobExecutorStudyIndexes.setIncludeIndexes(studyIndexes);
			rowJobExecutorStudyIndexes.setWriteFolder(studyBasedResultsFolder);

			//add the jobExecutors to the executor list
			rowJobExecutors.add(rowJobExecutorBiosIndexes);
			rowJobExecutors.add(rowJobExecutorStudyIndexes);

			double nanoTime = System.nanoTime();
			useExecutorsOnFile(	rowJobExecutors,
								correctedMatrixFn);
			double timeUsed = (System.nanoTime() - nanoTime) / 1000 / 1000 / 1000;
			System.out.println("timeused = " + timeUsed + " seconds");
			Thread.sleep(1000);
			closeExecutorFileWriters(rowJobExecutors);

			//OUTLIERCOUNTS Z-SCORES/////////////////////////////////////////////////////////////////////////////////////////////////////////
			//calculate outlierCounts in other parents
			//calculate outlierCounts in other patients
			String zScoreFn = studyBasedResultsFolder + biosBasedZscoreFn;
			HashMap<String, int[]> sampleToOtherParentIndexesStudy = getStudyIndexesPerSample(	zScoreFn,
																								stringFiltersStudy,
																								caseSampleNameToFamilySampleNames,
																								parentSampleNameToFamilyId);

			HashMap<String, int[]> sampleToOtherChildIndexesStudy = getStudyIndexesPerSample(	zScoreFn,
																								stringFiltersStudy,
																								caseSampleNameToFamilySampleNames,
																								patientSampleNameToFamilyId);

			ArrayList<RowJobExecutor> rowJobExecutors2 = new ArrayList<RowJobExecutor>();

			String cutoffPositiveCountFn = "parentOutlierCountsPositive.txt";
			RowAboveCutoffCounter rowAboveCutoffCounter = new RowAboveCutoffCounter();
			rowAboveCutoffCounter.setWriteFn(cutoffPositiveCountFn);
			rowAboveCutoffCounter.setCutoff(3);

			String belowNegativeCutoffCountFn = "parentOutlierCountsNegative.txt";
			RowBelowCutoffCounter rowBelowCutoffCounter = new RowBelowCutoffCounter();
			rowBelowCutoffCounter.setWriteFn(belowNegativeCutoffCountFn);
			rowBelowCutoffCounter.setCutoff(-3);

			for (String sampleName : sampleToOtherParentIndexesStudy.keySet())
			{
				String resultsFolder = this.writeFolder + sampleName + "/";
				FileUtils.makeDir(resultsFolder);
				RowJobExecutor rowJobExecutor = new RowJobExecutor();
				rowJobExecutor.addJob(rowAboveCutoffCounter);
				rowJobExecutor.addJob(rowBelowCutoffCounter);
				rowJobExecutor.setIncludeIndexes(sampleToOtherParentIndexesStudy.get(sampleName));
				rowJobExecutor.setWriteFolder(resultsFolder);
				rowJobExecutors2.add(rowJobExecutor);
			}

			String cutoffPositiveCountFnChildren = "childOutlierCountsPositive.txt";
			RowAboveCutoffCounter rowAboveCutoffCounterChildren = new RowAboveCutoffCounter();
			rowAboveCutoffCounterChildren.setWriteFn(cutoffPositiveCountFnChildren);
			rowAboveCutoffCounterChildren.setCutoff(3);

			String belowNegativeCutoffCountFnChildren = "childOutlierCountsNegative.txt";
			RowBelowCutoffCounter rowBelowCutoffCounterChildren = new RowBelowCutoffCounter();
			rowBelowCutoffCounterChildren.setWriteFn(belowNegativeCutoffCountFnChildren);
			rowBelowCutoffCounterChildren.setCutoff(-3);

			for (String sampleName : sampleToOtherChildIndexesStudy.keySet())
			{

				String resultsFolder = this.writeFolder + sampleName + "/";
				FileUtils.makeDir(resultsFolder);
				RowJobExecutor rowJobExecutor = new RowJobExecutor();
				rowJobExecutor.addJob(rowAboveCutoffCounterChildren);
				rowJobExecutor.addJob(rowBelowCutoffCounterChildren);
				rowJobExecutor.setIncludeIndexes(sampleToOtherChildIndexesStudy.get(sampleName));
				rowJobExecutor.setWriteFolder(resultsFolder);
				rowJobExecutors2.add(rowJobExecutor);
			}
			useExecutorsOnFile(	rowJobExecutors2,
								zScoreFn);
			closeExecutorFileWriters(rowJobExecutors2);
			
			//MERGE Z-SCORES FAMILY/////////////////////////////////////////////////////////////////////////////////////////////////////////
			//merge family files for z-scores based on BIOS
			ArrayList<RowJobExecutor> rowJobExecutors3 = new ArrayList<RowJobExecutor>();
			String zScoreFnFamily = "familyZscores.txt";

			RowColumnGetter columnGetter = new RowColumnGetter();
			columnGetter.setWriteFn(zScoreFnFamily);

			for (String sampleName : caseSampleNameToFamilySampleNames.keySet())
			{
				ArrayList<String> includeColnames = new ArrayList<String>();
				includeColnames.add(sampleName);
				includeColnames = addParentSamples(	caseSampleNameToFamilySampleNames,
													sampleName,
													includeColnames);
				//get family z_scores
				RowJobExecutor rowJobExecutor = new RowJobExecutor();
				rowJobExecutor.setIncludeColHeaders(includeColnames);
				rowJobExecutor.addJob(columnGetter);
				String resultsFolder = this.writeFolder + sampleName + "/";
				rowJobExecutor.setWriteFolder(resultsFolder);

				rowJobExecutors3.add(rowJobExecutor);
			}
			useExecutorsOnFile(	rowJobExecutors3,
								studyBasedResultsFolder + biosBasedZscoreFn);
			closeExecutorFileWriters(rowJobExecutors3);
			
			//RAWDATA/////////////////////////////////////////////////////////////////////////////////////////////////////////
			//calculate medians rawData BIOS
			String mediansFn = "medians_Bios.txt";
			RowMedianGetter rowMedianGetter = new RowMedianGetter();
			rowMedianGetter.setWriteFn(mediansFn);

			String maxRawCountBiosFn = "maxRawCountBios.txt";
			RowMaxGetter rowMaxGetter = new RowMaxGetter();
			rowMaxGetter.setWriteFn(maxRawCountBiosFn);

			String maxRawCountNonFamilyFn = "maxRawCountNonFamily.txt";
			RowMaxGetter rowMaxGetterNonFamily = new RowMaxGetter();
			rowMaxGetterNonFamily.setWriteFn(maxRawCountNonFamilyFn);

			String biosBasedResultsFolderRawData = this.writeFolder + "bios/rawData/";
			FileUtils.makeDir(biosBasedResultsFolderRawData);

			RowJobExecutor rowJobExecutorBiosIndexesRawData = new RowJobExecutor();
			rowJobExecutorBiosIndexesRawData.addJob(rowMedianGetter);
			rowJobExecutorBiosIndexesRawData.addJob(rowMaxGetter);
			rowJobExecutorBiosIndexesRawData.setIncludeIndexes(biosIndexes);
			rowJobExecutorBiosIndexesRawData.setWriteFolder(biosBasedResultsFolderRawData);

			String rawCountFn = "rawCount.txt";
			RowColumnGetter rawCountColumnGetter = new RowColumnGetter();
			rawCountColumnGetter.setWriteFn(rawCountFn);

			ArrayList<RowJobExecutor> rowJobExecutorsOnRawData = new ArrayList<RowJobExecutor>();
			for (String sampleName : caseSampleNameToFamilySampleNames.keySet())
			{
				//get raw counts
				RowJobExecutor rowJobExecutorGetRawCounts = new RowJobExecutor();
				ArrayList<String> includeColnamesCaseSample = new ArrayList<String>();
				includeColnamesCaseSample.add(sampleName);
				rowJobExecutorGetRawCounts.setIncludeColHeaders(includeColnamesCaseSample);
				rowJobExecutorGetRawCounts.addJob(rawCountColumnGetter);
				String resultsFolder = this.writeFolder + sampleName + "/";
				rowJobExecutorGetRawCounts.setWriteFolder(resultsFolder);

				rowJobExecutorsOnRawData.add(rowJobExecutorGetRawCounts);
			}

			//get max raw counts non family members
			for (String sampleName : sampleToOtherIndexesStudy.keySet())
			{
				int[] otherStudyIndexes = sampleToOtherIndexesStudy.get(sampleName);
				RowJobExecutor rowJobExecutorGetRawCounts = new RowJobExecutor();
				ArrayList<String> includeColnamesCaseSample = new ArrayList<String>();
				includeColnamesCaseSample.add(sampleName);
				rowJobExecutorGetRawCounts.setIncludeIndexes(otherStudyIndexes);

				rowJobExecutorGetRawCounts.addJob(rowMaxGetterNonFamily);
				String resultsFolder = this.writeFolder + sampleName + "/";
				rowJobExecutorGetRawCounts.setWriteFolder(resultsFolder);

				rowJobExecutorsOnRawData.add(rowJobExecutorGetRawCounts);
			}

			rowJobExecutorsOnRawData.add(rowJobExecutorBiosIndexesRawData);
			useExecutorsOnFile(	rowJobExecutorsOnRawData,
								this.rawMatrixFn);
			closeExecutorFileWriters(rowJobExecutorsOnRawData);

			///////////////////////////////////////////////////////////////////////////////////////////////////////////
			//MergeFiles and write Patrick files
			for (String caseSampleName : caseSampleNameToFamilySampleNames.keySet())
			{
				String familyName = getFamilyName(	caseSampleName,
													patientSampleNameToFamilyId,
													parentSampleNameToFamilyId);
				if (familyName == null)
					continue;
				String gnScoreFn = gnFolder + familyName + "_fixed.txt";
				HashMap<String, String> geneToGnScores = FileUtils.readStringStringHash(gnScoreFn,
																						false,
																						2,
																						1,
																						9);

				String vcfFn = vcfFolder + familyName + ".txt";
				HashMap<String, ArrayList<String>> ensgToPositions = FileUtils.readGavinVcfHash(vcfFn,
																								enstToEnsgFn);

				String sampleFolder = writeFolder + caseSampleName + "/";
				BufferedReader zScoreFamilyReader = FileUtils.createReader(sampleFolder + zScoreFnFamily);
				BufferedReader biosMedianReader = FileUtils.createReader(biosBasedResultsFolderRawData + mediansFn);

				BufferedReader otherParentsAboveCutoffCountReader = FileUtils.createReader(sampleFolder + cutoffPositiveCountFn);
				BufferedReader otherParentsBelowCutoffCountReader = FileUtils.createReader(sampleFolder + belowNegativeCutoffCountFn);

				BufferedReader caseAboveCutoffCountReader = FileUtils.createReader(sampleFolder + cutoffPositiveCountFnChildren);
				BufferedReader caseBelowCutoffCountReader = FileUtils.createReader(sampleFolder + belowNegativeCutoffCountFnChildren);

				BufferedReader rawCountReader = FileUtils.createReader(sampleFolder + rawCountFn);
				BufferedReader maxRawCountBiosReader = FileUtils.createReader(biosBasedResultsFolderRawData + maxRawCountBiosFn);
				BufferedReader maxRawCountNonFamilyReader = FileUtils.createReader(sampleFolder + maxRawCountNonFamilyFn);

				String mergeWriteFn = sampleFolder + "merged.txt";
				BufferedWriter mergedWriter = FileUtils.createWriter(mergeWriteFn);

				String patrickWriteFn = sampleFolder + caseSampleName + "_patrick.txt";
				BufferedWriter patrickWriter = FileUtils.createWriter(patrickWriteFn);

				String patrickFolder = writeFolder + "patrickFolder/";
				FileUtils.makeDir(patrickFolder);
				String patrickFolderWriteFn = patrickFolder + caseSampleName + ".txt";
				BufferedWriter patrickFolderWriter = FileUtils.createWriter(patrickFolderWriteFn);

				String rescueWriteFn = sampleFolder + caseSampleName + "_rescue.txt";
				BufferedWriter rescueWriter = FileUtils.createWriter(rescueWriteFn);

				String rescueFolder = writeFolder + "rescueFolder/";
				FileUtils.makeDir(rescueFolder);
				String rescueFolderWriteFn = rescueFolder + caseSampleName + ".txt";
				BufferedWriter rescueFolderWriter = FileUtils.createWriter(rescueFolderWriteFn);

				String line = null;
				boolean isHeaderLine = true;
				int l = 0;
				while ((line = zScoreFamilyReader.readLine()) != null)//assumes all files have the same rowNames in the same order
				{
					String[] rowName_zScoreFamily = line.split(	"\t",
																2);
					String rowName = rowName_zScoreFamily[0];
					String zScoreFamily = rowName_zScoreFamily[1];

					String biosMedian = biosMedianReader.readLine().split(	"\t",
																			2)[1];
					String otherParentsAboveCutoffCount = otherParentsAboveCutoffCountReader.readLine().split(	"\t",
																												2)[1];
					String caseAboveCutoffCount = caseAboveCutoffCountReader.readLine().split(	"\t",
																								2)[1];

					String otherParentsBelowCutoffCount = otherParentsBelowCutoffCountReader.readLine().split(	"\t",
																												2)[1];
					String caseBelowCutoffCount = caseBelowCutoffCountReader.readLine().split(	"\t",
																								2)[1];

					String rawCountStr = rawCountReader.readLine().split(	"\t",
																			2)[1];

					String maxRawCountBios = maxRawCountBiosReader.readLine().split("\t",
																					2)[1];
					String maxRawCountNonFamily = maxRawCountNonFamilyReader.readLine().split(	"\t",
																								2)[1];

					String mergedLine = getMergedLine(	rowName,
														zScoreFamily,
														biosMedian,
														otherParentsAboveCutoffCount,
														caseAboveCutoffCount,
														rawCountStr,
														maxRawCountBios,
														maxRawCountNonFamily);
					mergedWriter.write(mergedLine);

					writeCriteriaDependantLines(patrickWriter,
												patrickFolderWriter,
												rescueWriter,
												rescueFolderWriter,
												mergedLine,
												ensgNameToGeneSymbol,
												isHeaderLine,
												geneToGnScores,
												ensgToPositions,
												biosMedian,
												otherParentsAboveCutoffCount,
												caseAboveCutoffCount,
												zScoreFamily,
												rawCountStr,
												maxRawCountBios,
												maxRawCountNonFamily,
												otherParentsBelowCutoffCount,
												caseBelowCutoffCount);

					isHeaderLine = false;
					l++;
				}
				patrickWriter.close();
				patrickFolderWriter.close();
				rescueWriter.close();
				rescueFolderWriter.close();
				mergedWriter.close();

			}

		} catch (Exception e)
		{
			e.printStackTrace();
		}
	}

	private void writeRescueLine(	BufferedWriter rescueWriter,
									String mergedLine)
	{

	}

	private String getMergedLine(	String rowName,
									String zScoreFamily,
									String biosMedian,
									String biosCutoffCount,
									String caseCutoffCount,
									String rawCount,
									String maxRawCountBios,
									String maxRawCountNonFamily) throws IOException
	{
		StringBuilder stringBuilder = new StringBuilder();
		stringBuilder.append(rowName);
		stringBuilder.append("\t");
		stringBuilder.append(zScoreFamily);
		stringBuilder.append("\t");
		stringBuilder.append(biosMedian);
		stringBuilder.append("\t");
		stringBuilder.append(biosCutoffCount);
		stringBuilder.append("\t");
		stringBuilder.append(caseCutoffCount);
		stringBuilder.append("\t");
		stringBuilder.append(rawCount);
		stringBuilder.append("\t");
		stringBuilder.append(maxRawCountBios);
		stringBuilder.append("\t");
		stringBuilder.append(maxRawCountNonFamily);
		stringBuilder.append("\n");
		String writeLine = stringBuilder.toString();
		return writeLine;
	}

	private String getFamilyName(	String sampleName,
									HashMap<String, String> patientSampleNameToFamilyId,
									HashMap<String, String> parentSampleNameToFamilyId)
	{
		String familyName = patientSampleNameToFamilyId.get(sampleName);
		if (familyName == null)
		{
			familyName = parentSampleNameToFamilyId.get(sampleName);
			if (familyName != null)
			{
				p("Sample is parent: " + sampleName + ". Skipping");
				return null;
			}
		}
		p("familyName =\t" + familyName);
		if (familyName == null)
		{
			p("Sample to family name missing for sample: " + sampleName + ". Skipping");
			return null;
		}
		return familyName;
	}

	private void writeCriteriaDependantLines(	BufferedWriter patrickWriter,
												BufferedWriter patrickFolderWriter,
												BufferedWriter rescueWriter,
												BufferedWriter rescueFolderWriter,
												String line,
												HashMap<String, String> ensgNameToGeneSymbol,
												boolean isHeaderLine,
												HashMap<String, String> geneToGnScores,
												HashMap<String, ArrayList<String>> ensgToPositions,
												String biosMedian,
												String otherparentsCutoffCount,
												String caseCutoffCount,
												String zScoreFamily,
												String rawCountStr,
												String maxRawCountBios,
												String maxRawCountNonFamily,
												String otherParentsBelowCutoffCount,
												String caseBelowCutoffCount) throws IOException
	{
		if (!isHeaderLine)
		{

			double median = Double.parseDouble(biosMedian);
			double outlierCountParents = Double.parseDouble(otherparentsCutoffCount);
			double outlierCountCases = Double.parseDouble(caseCutoffCount);
			
			double negativeOutlierCountParents = Double.parseDouble(otherParentsBelowCutoffCount);
			double negativeOutlierCountCases = Double.parseDouble(caseBelowCutoffCount);
			
			String patientZscoreAsString = zScoreFamily.split("\t")[0];
			double patientZscore = Double.parseDouble(patientZscoreAsString);
			String parent1ZscoreAsString = zScoreFamily.split("\t")[1];
			String parent2ZscoreAsString = zScoreFamily.split("\t")[2];
			double rawCount = Double.parseDouble(rawCountStr);
			double maxBiosRawCount = Double.parseDouble(maxRawCountBios);
			double maxStudyRawCount = Double.parseDouble(maxRawCountNonFamily);

			//get the variables that are needed to determine if the line is passing the criteria
			//I should get the indexes from the columnheaders but cant be bothered at the moment...
			String[] ensembl_spliceName = line.split(	"\t",
														2)[0].split("__");
			String ensemblName = ensembl_spliceName[0];
			String spliceName = ensembl_spliceName[1];
			if (ensemblName.equals("null"))
			{
				ensemblName = "unannotated";
			}

			String[] valuesAsString = line.split(	"\t",
													2)[1].split("\t");

			String geneSymbol = getFromHash(ensgNameToGeneSymbol,
											ensemblName);
			String geneNetworkScore = getFromHash(	geneToGnScores,
													ensemblName);

			patrickCriteriaPassWrite(	outlierCountParents,
										outlierCountCases,
										patientZscore,
										ensgNameToGeneSymbol,
										ensemblName,
										geneToGnScores,
										patrickWriter,
										patrickFolderWriter,
										ensgToPositions,
										spliceName,
										median,
										patientZscoreAsString,
										parent1ZscoreAsString,
										parent2ZscoreAsString,
										geneSymbol,
										geneNetworkScore,
										negativeOutlierCountParents,
										negativeOutlierCountCases);

			rescueCriteriaPassWrite(outlierCountParents,
									outlierCountCases,
									patientZscore,
									ensgNameToGeneSymbol,
									ensemblName,
									geneToGnScores,
									rescueWriter,
									rescueFolderWriter,
									ensgToPositions,
									spliceName,
									median,
									patientZscoreAsString,
									parent1ZscoreAsString,
									parent2ZscoreAsString,
									geneSymbol,
									geneNetworkScore,
									rawCount,
									maxBiosRawCount,
									maxStudyRawCount);

		}
		else
		{
			String patrickHeader = "Gene	ENSG	geneNetwork Z-score	Case expression	Father expression	Mother expression	Variant chr	Variant pos	Variant ID	Ref Allele	Alt Allele	Case genotype	Exac MAF\n";
			patrickWriter.write(patrickHeader);
			patrickFolderWriter.write(patrickHeader);
			String rescueHeader = "Gene	ENSG	geneNetwork Z-score	Case expression	Father expression	Mother expression	Variant chr	Variant pos	Variant ID	Ref Allele	Alt Allele	Case genotype	Exac MAF\n";
			rescueWriter.write(rescueHeader);
			rescueFolderWriter.write(rescueHeader);
		}
	}

	private void rescueCriteriaPassWrite(	double outlierCountParents,
											double outlierCountCases,
											double patientZscore,
											HashMap<String, String> ensgNameToGeneSymbol,
											String ensemblName,
											HashMap<String, String> geneToGnScores,
											BufferedWriter rescueWriter,
											BufferedWriter rescueFolderWriter,
											HashMap<String, ArrayList<String>> ensgToPositions,
											String spliceName,
											double median,
											String patientZscoreAsString,
											String parent1ZscoreAsString,
											String parent2ZscoreAsString,
											String geneSymbol,
											String geneNetworkScore,
											double rawCount,
											double maxBiosRawCount,
											double maxStudyRawCount) throws IOException
	{
		boolean isMuchLargerThanBios = rawCount > maxBiosRawCount * 5 || maxBiosRawCount == 0;
		boolean isMuchLargerThanStudy = rawCount > maxStudyRawCount * 2;
		if (median <= medianCutoff && isMuchLargerThanBios && isMuchLargerThanStudy && Math.abs(patientZscore) > zScoreCutoff)
		{
			String rescueLine = geneSymbol.concat("\t").concat(ensemblName).concat("__").concat(spliceName).concat("\t").concat(geneNetworkScore).concat("\t").concat(patientZscoreAsString).concat("\t").concat(parent1ZscoreAsString).concat("\t").concat(parent2ZscoreAsString);

			addVariantsAndWriteLine(rescueLine,
									rescueWriter,
									rescueFolderWriter,
									ensgToPositions,
									ensemblName);

		}
	}

	private void patrickCriteriaPassWrite(	double outlierCountParents,
											double outlierCountCases,
											double patientZscore,
											HashMap<String, String> ensgNameToGeneSymbol,
											String ensemblName,
											HashMap<String, String> geneToGnScores,
											BufferedWriter patrickWriter,
											BufferedWriter patrickFolderWriter,
											HashMap<String, ArrayList<String>> ensgToPositions,
											String spliceName,
											double median,
											String patientZscoreAsString,
											String parent1ZscoreAsString,
											String parent2ZscoreAsString,
											String geneSymbol,
											String geneNetworkScore,
											double negativeOutlierCountParents,
											double negativeOutlierCountCases) throws IOException
	{
		boolean isMedianPass = median > medianCutoff;
		boolean isPositivePass = outlierCountParents < maxOutlierCount && outlierCountCases < maxOutlierCount && patientZscore > zScoreCutoff;
		boolean isNegativepass = negativeOutlierCountParents < maxOutlierCount && negativeOutlierCountCases < maxOutlierCount && patientZscore < -zScoreCutoff;

		if (isMedianPass && (isPositivePass || isNegativepass))
		{
			String patrickLine = geneSymbol.concat("\t").concat(ensemblName).concat("__").concat(spliceName).concat("\t").concat(geneNetworkScore).concat("\t").concat(patientZscoreAsString).concat("\t").concat(parent1ZscoreAsString).concat("\t").concat(parent2ZscoreAsString);

			addVariantsAndWriteLine(patrickLine,
									patrickWriter,
									patrickFolderWriter,
									ensgToPositions,
									ensemblName);

		}
	}

	private void addVariantsAndWriteLine(	String patrickLine,
											BufferedWriter patrickWriter,
											BufferedWriter patrickFolderWriter,
											HashMap<String, ArrayList<String>> ensgToPositions,
											String ensemblName) throws IOException
	{
		String writeLine = patrickLine;

		ArrayList<String> positions = null;
		if (ensgToPositions != null)
		{
			positions = ensgToPositions.get(ensemblName);
		}

		if (positions != null)
			for (String position : positions)
			{
				writeLine = patrickLine;
				String[] eles = position.split("\t");
				writeLine = writeLine.concat("\t");
				writeLine = writeLine.concat(eles[0]);//chromosome
				writeLine = writeLine.concat("\t");
				writeLine = writeLine.concat(eles[1]);//position
				writeLine = writeLine.concat("\t");
				writeLine = writeLine.concat(eles[2]);//rs number
				writeLine = writeLine.concat("\t");
				writeLine = writeLine.concat("\t");
				writeLine = writeLine.concat(eles[4]);//alt allele

				//case genotype
				String genotype = eles[5].split(":")[0];
				writeLine = writeLine.concat("\t");
				writeLine = writeLine.concat(genotype);

				//exacMaf
				String exacMaf = eles[7];
				writeLine = writeLine.concat("\t");
				writeLine = writeLine.concat(exacMaf);
				writeLine = writeLine.concat("\n");
				patrickWriter.write(writeLine);
				if (patrickFolderWriter != null)
					patrickFolderWriter.write(writeLine);
			}
		else
		{
			for (int i = 0; i < 7; i++)
			{
				writeLine = writeLine.concat("\t");
				writeLine = writeLine.concat("");
			}
			writeLine = writeLine.concat("\n");
			patrickWriter.write(writeLine);
			if (patrickFolderWriter != null)
				patrickFolderWriter.write(writeLine);
		}
	}

	private ArrayList<String> addParentSamples(	HashMap<String, HashMap<String, Integer>> sampleNameToFamilySampleNames,
												String sampleName,
												ArrayList<String> includeColnames)
	{
		HashMap<String, Integer> includeSamples = sampleNameToFamilySampleNames.get(sampleName);
		includeSamples.remove(sampleName);
		for (String includeSample : includeSamples.keySet())
		{
			includeColnames.add(includeSample);
		}
		return includeColnames;
	}

	private void closeExecutorFileWriters(ArrayList<RowJobExecutor> rowJobExecutors) throws IOException
	{
		for (RowJobExecutor rowJobExecutor : rowJobExecutors)
		{
			rowJobExecutor.closeWriters();
		}
	}

	private void useExecutorsOnFile(ArrayList<RowJobExecutor> rowJobExecutors,
									String matrixFn) throws FileNotFoundException, IOException, InterruptedException
	{
		BufferedReader sampleReader = FileUtils.createReader(matrixFn);
		String header = sampleReader.readLine();//get rid of header

		executeHeaderJobs(	rowJobExecutors,
							header);

		//		int queueSize = nThreads * 2;
		//		ExecutorService executorService = new ThreadPoolExecutor(	nThreads,
		//																	nThreads,
		//																	5000L,
		//																	TimeUnit.MILLISECONDS,
		//																	new ArrayBlockingQueue<Runnable>(	10,
		//																										true),
		//																	new ThreadPoolExecutor.CallerRunsPolicy());
		//		Executors.newFixedThreadPool(nThreads);

		String line = null;
		int lineNumber = 0;
		while ((line = sampleReader.readLine()) != null)
		{
			//			Runnable worker = new RowJobExecutorThread(	line,
			//														lineNumber,
			//														rowJobExecutors);
			//			lineNumber++;
			//			executorService.execute(worker);

			String rowName = line.split("\t",
										2)[0];
			double[] values =null;

			try
			{
			String[] valuesString = line.split(	"\t",
												2)[1].split("\t");
			

				values = FileUtils.convertToDoubleArray(valuesString);
			}catch(Exception e)
			{
				e.printStackTrace();
				System.out.println(line);
				System.out.println(matrixFn);
			}

			setRowJobExecutorVariables(	rowJobExecutors,
										rowName,
										values,
										lineNumber);

			executeJobs(rowJobExecutors);

		}
		//		executorService.shutdown();
		//		while (!executorService.isTerminated())
		//		{
		//			Thread.sleep(1);
		//		}
		//		;
	}

	private void setRowJobExecutorVariables(ArrayList<RowJobExecutor> rowJobExecutors,
											String rowName,
											double[] values,
											int lineNumber)
	{
		for (RowJobExecutor rowJobExecutor : rowJobExecutors)
		{
			rowJobExecutor.setRowName(rowName);

			setValuesUsingIncludeIndexes(	values,
											rowJobExecutor);
		}
	}

	public void setValuesUsingIncludeIndexes(	double[] values,
												RowJobExecutor rowJobExecutor)
	{
		rowJobExecutor.initiateStorageVariables();
		int[] includeIndexes = rowJobExecutor.getIncludeIndexes();
		double[] includeValues = new double[includeIndexes.length];
		for (int i = 0; i < includeIndexes.length; i++)
		{
			int e = includeIndexes[i];
			includeValues[i] = values[e];
		}
		rowJobExecutor.setValuesUsingIncludeIndexes(values);
	}

	private void executeJobs(ArrayList<RowJobExecutor> rowJobExecutors)
	{
		for (RowJobExecutor rowJobExecutor : rowJobExecutors)
		{
			rowJobExecutor.execute(0);
		}
	}

	private void executeHeaderJobs(	ArrayList<RowJobExecutor> rowJobExecutors,
									String header)
	{
		for (RowJobExecutor rowJobExecutor : rowJobExecutors)
		{
			rowJobExecutor.executeHeaderJob(header);
		}
	}

	private void createWriteFolder()
	{
		if (writeFolder == null)
			writeFolder = new File(correctedMatrixFn).getParent() + "/R.SplicingAnalysisPipeline2/";
		writeFolder = FileUtils.makeFolderNameEndWithSlash(writeFolder);
		new File(writeFolder).mkdir();
	}

	private int[] getIndexes(	String correctedMatrixFn2,
								StringFilters stringFilters) throws FileNotFoundException, IOException
	{
		return getIndexes(	correctedMatrixFn2,
							stringFilters,
							null,
							null);
	}

	private int[] getIndexes(	String matrixFn,
								StringFilters stringFilters,
								HashMap<String, Integer> excludeSamples,
								HashMap<String, String> includeSamples) throws FileNotFoundException, IOException
	{
		BufferedReader reader = FileUtils.createReader(matrixFn);
		String header = reader.readLine();
		String[] colNames = header.split(	"\t",
											2)[1].split("\t");//first split the first column off
		boolean[] include = new boolean[colNames.length];
		int n = 0;
		for (int c = 0; c < colNames.length; c++)
		{
			include[c] = stringFilters.check(colNames[c]);
			//exclude family samples			
			if (excludeSamples != null && excludeSamples.containsKey(colNames[c]) || (includeSamples != null && !includeSamples.containsKey(colNames[c])))
			{
				include[c] = false;
			}
			if (include[c])
				n++;
		}
		int[] biosIndexes = new int[n];
		int arrayIndex = 0;
		for (int b = 0; b < include.length; b++)
		{
			if (!include[b])
				continue;

			biosIndexes[arrayIndex] = b;
			arrayIndex++;
		}

		reader.close();
		return biosIndexes;
	}

	private HashMap<String, int[]> getStudyIndexesPerSample(String correctedMatrixFn,
															StringFilters stringFilters,
															HashMap<String, HashMap<String, Integer>> sampleNameToFamilySampleNames) throws FileNotFoundException, IOException
	{
		return getStudyIndexesPerSample(correctedMatrixFn,
										stringFilters,
										sampleNameToFamilySampleNames,
										null);
	}

	private HashMap<String, int[]> getStudyIndexesPerSample(String matrixFn,
															StringFilters stringFilters,
															HashMap<String, HashMap<String, Integer>> sampleNameToFamilySampleNames,
															HashMap<String, String> includeSamples) throws FileNotFoundException, IOException
	{

		BufferedReader reader = FileUtils.createReader(matrixFn);
		String header = reader.readLine();
		String[] colNames = header.split(	"\t",
											2)[1].split("\t");//first split the first column off
		reader.close();
		HashMap<String, int[]> studyIndexesPerSample = new HashMap<String, int[]>();

		for (int c = 0; c < colNames.length; c++)
		{

			if (stringFilters.check(colNames[c]))
			{
				String sampleName = colNames[c];
				HashMap<String, Integer> familySampleNames = sampleNameToFamilySampleNames.get(sampleName);

				if (familySampleNames != null)
				{
					int[] sampleIncludeIndexes = getIndexes(matrixFn,
															stringFilters,
															familySampleNames,
															includeSamples);

					studyIndexesPerSample.put(	sampleName,
												sampleIncludeIndexes);
				}
				else//is probably not a child, but a parent
					;//p("Warning, sampleName: " + sampleName + " not found in sampleNameToFamilySampleNames file:\t" + this.sampleNameToFamilySampleNamesFn);

			}
		}

		return studyIndexesPerSample;
	}

	private String getFromHash(	HashMap<String, String> geneToMedianBios,
								String geneName) throws IOException
	{
		String median = null;
		if (geneToMedianBios != null)
		{
			median = geneToMedianBios.get(geneName);
		}
		if (median == null)
		{
			//System.out.println("Warning key "+geneName+" not found, using -1");
			median = "";
		}
		return median;
	}
}
