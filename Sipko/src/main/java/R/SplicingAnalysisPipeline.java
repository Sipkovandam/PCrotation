package R;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;

import org.apache.commons.math3.stat.descriptive.rank.Median;

import MatrixScripts.MergeFiles;
import MatrixScripts.Row;
import MatrixScripts.RowStatisticsGetter;
import MatrixScripts.ZscoreCalculator;
import Tools.FileSearcher;
import Tools.FileUtils;
import Tools.Script;
import Tools.StringFilters;
import Tools.StringFilters.ForbiddenStringChecker;

public class SplicingAnalysisPipeline extends Script<SplicingAnalysisPipeline>
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
			HashMap<String, HashMap<String, Integer>> sampleNameToFamilySampleNames = FileUtils.readStringMultiStringHash(	sampleNameToFamilySampleNamesFn,
																															0,
																															1,
																															false);

			HashMap<String, Double> maxBiosRawCountTracker = new HashMap<String, Double>();
			HashMap<String, String> spliceToMedianBios = new HashMap<String, String>();

			HashMap<String, String> ensgNameToGeneSymbol = FileUtils.readStringStringHash(	ensgNameToGeneSymbolFn,
																							false,
																							1,
																							1,
																							0);

			ArrayList<String> foldersToCombine = new ArrayList<String>();

			createWriteFolder();

			//make z-score tables for corrected files
			foldersToCombine.addAll(createZscoresCorrectedFiles(patientSampleNameToFamilyId,
																parentSampleNameToFamilyId,
																sampleNameToFamilySampleNames,
																writeFolder));

			//make z-score tables for uncorrected files
			foldersToCombine.addAll(createZscoresUncorrectedFiles(	patientSampleNameToFamilyId,
																	parentSampleNameToFamilyId,
																	sampleNameToFamilySampleNames,
																	writeFolder,
																	maxBiosRawCountTracker,
																	spliceToMedianBios));

			//combine all together

			File[] files = new File(foldersToCombine.get(0)).listFiles();
			String writeFolderCombined = this.writeFolder + "combinedResults/";
			String writeFolderSummary = this.writeFolder + "combinedResultsSummary/";
			String writeFolderPatrick = this.writeFolder + "combinedResultsPatrick/";

			new File(writeFolderCombined).mkdir();
			new File(writeFolderSummary).mkdir();
			new File(writeFolderPatrick).mkdir();

			for (File file : files)
			{
				String sampleFile = file.getName();
				String sampleName = FileUtils.removeExtention(sampleFile);
				p("sampleName =\t" + sampleFile);
				//get the family name
				String familyName = patientSampleNameToFamilyId.get(sampleName);
				if (familyName == null)
				{
					familyName = parentSampleNameToFamilyId.get(sampleName);
					if (familyName != null)
					{
						p("Sample is parent: " + sampleName + ". Skipping");
						continue;
					}
				}
				p("familyName =\t" + familyName);
				if (familyName == null)
				{
					p("Sample to family name missing for sample: " + sampleName + ". Skipping");
					continue;
				}
				//get the gn Scores
				String gnScoreFn = gnFolder + familyName + "_fixed.txt";
				HashMap<String, String> geneToGnScores = FileUtils.readStringStringHash(gnScoreFn,
																						false,
																						2,
																						1,
																						9);
				//get the mutations
				String vcfFn = vcfFolder + familyName + ".txt";
				HashMap<String, ArrayList<String>> ensgToPositions = FileUtils.readGavinVcfHash(vcfFn,
																								enstToEnsgFn);

				BufferedWriter combinedWriter = FileUtils.createWriter(writeFolderCombined + sampleFile);
				BufferedWriter summaryWriter = FileUtils.createWriter(writeFolderSummary + FileUtils.addBeforeExtention(sampleFile,
																														"_summary"));
				BufferedWriter patrickWriter = FileUtils.createWriter(writeFolderPatrick + FileUtils.addBeforeExtention(sampleFile,
																														"_Patrick"));
				ArrayList<BufferedReader> sampleReaders = createReaders(foldersToCombine,
																		sampleFile);

				HashMap<String, Integer> parentIds = sampleNameToFamilySampleNames.get(sampleName);

				//Create z-score readers for the parents of the patient (if this sample is a patient)
				ArrayList<BufferedReader> parentZscoreReaders = new ArrayList<BufferedReader>();
				if (parentIds != null)
					for (String parentId : parentIds.keySet())
					{
						if (parentId.compareTo(sampleName) == 0)
							continue;
						String parentZscoreFn = foldersToCombine.get(0) + parentId + ".txt";
						p("parentZscoreFn" + parentZscoreFn);
						BufferedReader parentReader = FileUtils.createReader(parentZscoreFn);
						//parentReader.readLine();//remove/skip header
						parentZscoreReaders.add(parentReader);
					}

				String line = null;
				boolean firstLine = true;
				String[] parentBiosZscores = new String[parentZscoreReaders.size()];

				//go through the 8 files line by line simultaneously
				while ((line = sampleReaders.get(0).readLine()) != null)
				{
					String[] lineEles = line.split("\t");
					String newLine = "";
					String geneName = lineEles[0];
					String geneWithSpliceName = geneName;
					String spliceName = geneName;

					if (geneName.contains("__"))
					{
						String[] geneAndSplice = geneName.split("__");
						geneName = geneAndSplice[0];
						spliceName = geneAndSplice[1];
						if (geneName.matches("null") || geneName.matches(""))
						{
							geneName = "unannotated";
							geneWithSpliceName = geneWithSpliceName.replaceAll(	".*__",
																				geneName + "__");
						}
					}

					String biosZscore = lineEles[1];
					newLine = newLine.concat(line);

					//get one column from each file
					for (int c = 1; c < sampleReaders.size(); c++)
					{
						newLine = newLine.concat("\t");
						newLine = newLine.concat(sampleReaders.get(c).readLine().split(	"\t",
																						2)[1]);
					}

					for (int r = 0; r < parentZscoreReaders.size(); r++)
					{
						//get the z-score (2nd column)
						String parentLine = parentZscoreReaders.get(r).readLine();
						String zScoreParent = parentLine.split("\t")[1];
						parentBiosZscores[r] = zScoreParent;
					}

					newLine = addFromHash(	spliceToMedianBios,
											spliceName,
											"Median_raw_expression_Bios",
											firstLine,
											newLine);

					//write GeneNetwork scores
					newLine = addFromHash(	geneToGnScores,
											geneName,
											"GeneNetworkScore",
											firstLine,
											newLine);

					//write medians
					double median = 0;
					if (firstLine == false)
						median = Double.parseDouble(getFromHash(spliceToMedianBios,
																spliceName));

					ArrayList<String> positions = null;
					if (ensgToPositions != null)
					{
						positions = ensgToPositions.get(geneName);
					}
					double maxBios = maxBiosRawCountTracker.get(spliceName);
					if (positions != null)
						for (String position : positions)
						{
							String writeLine = newLine;
							writeLine = writeLine.concat("\t");

							writeLine = writeLine.concat(position);
							writeLine = writeLine.concat("\n");
							combinedWriter.write(writeLine);

							writeSummaryFiles(	writeLine,
												patientSampleNameToFamilyId,
												sampleName,
												firstLine,
												summaryWriter,
												patrickWriter,
												parentZscoreReaders,
												biosZscore,
												ensgNameToGeneSymbol,
												geneName,
												geneToGnScores,
												parentBiosZscores,
												position,
												geneWithSpliceName,
												median,
												maxBios);
						}
					else
					{
						String writeLine = newLine;
						writeLine = writeLine.concat("\t");

						writeLine = writeLine.concat("");
						writeLine = writeLine.concat("\n");
						combinedWriter.write(writeLine);

						writeSummaryFiles(	writeLine,
											patientSampleNameToFamilyId,
											sampleName,
											firstLine,
											summaryWriter,
											patrickWriter,
											parentZscoreReaders,
											biosZscore,
											ensgNameToGeneSymbol,
											geneName,
											geneToGnScores,
											parentBiosZscores,
											null,
											geneWithSpliceName,
											median,
											maxBios);
					}
					firstLine = false;
				}
				summaryWriter.close();
				combinedWriter.close();
				patrickWriter.close();
				closeStreams(sampleReaders);
			}
		} catch (Exception e)
		{
			e.printStackTrace();
		}

	}

	private void createWriteFolder()
	{
		if (writeFolder == null)
			writeFolder = new File(correctedMatrixFn).getParent() + "/R.SplicingAnalysisPipeline/";
		writeFolder = FileUtils.makeFolderNameEndWithSlash(writeFolder);
		new File(writeFolder).mkdir();
	}

	private void writeSummaryFiles(	String newLine,
									HashMap<String, String> patientSampleNameToFamilyId,
									String sampleName,
									boolean firstLine,
									BufferedWriter summaryWriter,
									BufferedWriter patrickWriter,
									ArrayList<BufferedReader> parentZscoreReaders,
									String biosZscore,
									HashMap<String, String> ensgNameToGeneSymbol,
									String geneName,
									HashMap<String, String> geneToGnScores,
									String[] parentBiosZscores,
									String position,
									String spliceName,
									double median,
									double maxBios) throws IOException
	{
		if (patientSampleNameToFamilyId.containsKey(sampleName))
		{
			if (firstLine)
			{
				summaryWriter.write(newLine);
				writePatrickHeader(	patrickWriter,
									parentZscoreReaders);
			}
			else
			{
				boolean isFailCriteria1 = median < medianCutoff;
				boolean isFailCriteria2 = median >= medianCutoff && median > 5 * maxBios;

				if (isFailCriteria1 || isFailCriteria2)
					return;
				String[] newEles = newLine.split("\t");
				Double zScore = Double.parseDouble(biosZscore);
				if (Double.isNaN(zScore) || Math.abs(zScore) > summaryCutoff)
				{
					if (!isMultiOutlier(newEles))
					{
						summaryWriter.write(newLine);

						writePatrickLine(	ensgNameToGeneSymbol,
											geneName,
											geneToGnScores,
											firstLine,
											newLine,
											biosZscore,
											patrickWriter,
											parentBiosZscores,
											position,
											spliceName);
					}

				}
			}
		}

	}

	private void writePatrickLine(	HashMap<String, String> ensgNameToGeneSymbol,
									String geneName,
									HashMap<String, String> geneToGnScores,
									boolean firstLine,
									String newLine,
									String biosZscore,
									BufferedWriter patrickWriter,
									String[] parentBiosZscores,
									String position,
									String spliceName) throws IOException
	{

		String geneSymbol = ensgNameToGeneSymbol.get(geneName);
		if (geneSymbol == null)
			geneSymbol = "";
		String patrickLine = geneSymbol.concat("\t").concat(spliceName);
		patrickLine = addFromHash(	geneToGnScores,
									geneName,
									"GeneNetworkScore",
									firstLine,
									patrickLine);

		//add child z-score
		patrickLine = patrickLine.concat("\t");
		patrickLine = patrickLine.concat(biosZscore);
		//add parent z-scores
		for (String parentBiosZscore : parentBiosZscores)
		{
			patrickLine = patrickLine.concat("\t");
			patrickLine = patrickLine.concat(parentBiosZscore);
		}

		if (position != null)
		{
			String[] eles = position.split("\t");
			patrickLine = patrickLine.concat("\t");
			patrickLine = patrickLine.concat(eles[0]);//chromosome
			patrickLine = patrickLine.concat("\t");
			patrickLine = patrickLine.concat(eles[1]);//position
			patrickLine = patrickLine.concat("\t");
			patrickLine = patrickLine.concat(eles[2]);//rs number
			patrickLine = patrickLine.concat("\t");
			patrickLine = patrickLine.concat(eles[3]);//ref allele
			patrickLine = patrickLine.concat("\t");
			patrickLine = patrickLine.concat(eles[4]);//alt allele

			//case genotype
			String genotype = eles[5].split(":")[0];
			patrickLine = patrickLine.concat("\t");
			patrickLine = patrickLine.concat(genotype);

			//exacMaf
			String exacMaf = eles[7];
			patrickLine = patrickLine.concat("\t");
			patrickLine = patrickLine.concat(exacMaf);
		}
		else
		{
			for (int i = 0; i < 7; i++)
			{
				patrickLine = patrickLine.concat("\t");
				patrickLine = patrickLine.concat("");
			}
		}
		patrickLine = patrickLine.concat("\n");
		patrickWriter.write(patrickLine);
	}

	private void writePatrickHeader(BufferedWriter patrickWriter,
									ArrayList<BufferedReader> parentZscoreReaders) throws IOException
	{
		//"Gene	ENSG	GeneNetwork Z-score	Case expression	Father expression	Mother expression	Variant chr	Variant pos	Variant ID	Ref Allele	Alt Allele	Case genotype	Exac MAF	CGD condition	CGD inheritance	CGD age group	CGD manifestation	Predicted MOI	Posterior AD	Posterior AR\n"
		patrickWriter.write("Gene	ENSG	GeneNetwork Z-score	Case expression	Father expression	Mother expression	Variant chr	Variant pos	Variant ID	Ref Allele	Alt Allele	Case genotype	Exac MAF\n");
	}

	private ArrayList<String> createZscoresUncorrectedFiles(HashMap<String, String> patientSampleNameToFamilyId,
															HashMap<String, String> parentSampleNameToFamilyId,
															HashMap<String, HashMap<String, Integer>> sampleNameToFamilySampleNames,
															String rootFolder,
															HashMap<String, Double> maxBiosRawCount,
															HashMap<String, String> biosMedian) throws FileNotFoundException, IOException
	{
		String[] descriptionsUncorrected = new String[] { "Bios", "Study" };
		String writeFolder = rootFolder + FileUtils.removeExtention(new File(uncorrectedMatrixFn).getName());
		ArrayList<String> foldersUncorrected = new ArrayList<String>();
		foldersUncorrected.add(writeFolder + "_" + descriptionsUncorrected[0] + "/");
		foldersUncorrected.add(writeFolder + "_" + descriptionsUncorrected[1] + "/");

		String childOutlierFolderUnCorrected = writeFolder + "childOutliercounts/" + descriptionsUncorrected[0] + "/";
		String parentOutlierFolderUnCorrected = writeFolder + "parentOutliercounts/" + descriptionsUncorrected[0] + "/";

		String maxReadCountOtherCasesFolder = writeFolder + "maxReadCountOtherCases/" + descriptionsUncorrected[0] + "/";
		String maxReadCountOtherParentFolder = writeFolder + "maxReadCountOtherParent/" + descriptionsUncorrected[0] + "/";

		executeSteps(	uncorrectedMatrixFn,
						this.ourStudyStrings,
						descriptionsUncorrected,
						foldersUncorrected,
						sampleNameToFamilySampleNames,
						patientSampleNameToFamilyId,
						parentSampleNameToFamilyId,
						writeFolder,
						parentOutlierFolderUnCorrected,
						childOutlierFolderUnCorrected,
						maxReadCountOtherCasesFolder,
						maxReadCountOtherParentFolder,
						maxBiosRawCount,
						biosMedian);

		foldersUncorrected.add(maxReadCountOtherCasesFolder);
		foldersUncorrected.add(childOutlierFolderUnCorrected);

		foldersUncorrected.add(parentOutlierFolderUnCorrected);
		foldersUncorrected.add(maxReadCountOtherParentFolder);

		return foldersUncorrected;
	}

	private ArrayList<String> createZscoresCorrectedFiles(	HashMap<String, String> patientSampleNameToFamilyId,
															HashMap<String, String> parentSampleNameToFamilyId,
															HashMap<String, HashMap<String, Integer>> sampleNameToFamilySampleNames,
															String rootFolder) throws FileNotFoundException, IOException
	{
		String[] descriptionsCorrected = new String[] { "Bios", "Study" };
		ArrayList<String> foldersCorrected = new ArrayList<String>();
		String writeFolder = rootFolder + FileUtils.removeExtention(new File(correctedMatrixFn).getName());
		foldersCorrected.add(writeFolder + "_" + descriptionsCorrected[0] + "/");
		foldersCorrected.add(writeFolder + "_" + descriptionsCorrected[1] + "/");

		String parentOutlierFolderCorrected = writeFolder + "parentOutliercounts/" + descriptionsCorrected[0] + "/";
		String childOutlierFolderCorrected = writeFolder + "childOutliercounts/" + descriptionsCorrected[0] + "/";
		ArrayList<String> folders = new ArrayList<String>();
//		executeSteps(	correctedMatrixFn,
//						this.ourStudyStrings,
//						descriptionsCorrected,
//						foldersCorrected,
//						sampleNameToFamilySampleNames,
//						patientSampleNameToFamilyId,
//						parentSampleNameToFamilyId,
//						writeFolder,
//						parentOutlierFolderCorrected,
//						childOutlierFolderCorrected,
//						null,
//						null,
//						null,
//						null);

		foldersCorrected.add(parentOutlierFolderCorrected);
		foldersCorrected.add(childOutlierFolderCorrected);
		return foldersCorrected;
	}

	private boolean isMultiOutlier(String[] newEles)
	{
		int otherCaseOutlierCount = Integer.parseInt(newEles[3]);
		int otherParentOutlierCount = Integer.parseInt(newEles[4]);
		if (otherCaseOutlierCount > this.maxOutlierCount || otherParentOutlierCount > this.maxOutlierCount)
			return true;
		return false;
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
			median = "NaN";
		}
		return median;
	}

	private String addFromHash(	HashMap<String, String> geneToMedianBios,
								String geneName,
								String colHeader,
								boolean firstLine,
								String newLine) throws IOException
	{
		if (firstLine == true)
		{
			newLine = newLine.concat("\t");
			newLine = newLine.concat(colHeader);
		}
		else
		{
			String median = getFromHash(geneToMedianBios,
										geneName);
			newLine = newLine.concat("\t");
			newLine = newLine.concat(median);
		}
		return newLine;
	}

	private void calculateAndWriteMediansFromBios(	String matrixFn,
													String mediansWriteFn) throws FileNotFoundException, IOException
	{
		StringFilters stringFiltersBios = new StringFilters();
		stringFiltersBios.addFilter(stringFiltersBios.new ForbiddenStringChecker(ourStudyStrings));
		int[] biosIndexes = getIndexes(	matrixFn,
										stringFiltersBios);

		makeAndWriteMedians(matrixFn,
							biosIndexes,
							mediansWriteFn);
	}

	private void makeAndWriteMedians(	String matrixFn,
										int[] biosIndexes,
										String mediansWriteFn) throws FileNotFoundException, IOException
	{
		BufferedReader sampleReader = FileUtils.createReader(matrixFn);
		BufferedWriter mediansWriter = FileUtils.createWriter(mediansWriteFn);

		String header = sampleReader.readLine();
		mediansWriter.write(header + "\n");
		mediansWriter.write("RowName\tmedianBiosReads\n");

		String line = null;
		while ((line = sampleReader.readLine()) != null)
		{
			//read the row
			Row row = Row.readRow(line);
			double median = getMedian(row.getValues());
			//write these to file
			mediansWriter.write(row.getRowName().concat("\t").concat(Double.toString(median).concat("\n")));
		}
		mediansWriter.close();
	}

	private ArrayList<BufferedReader> createReaders(ArrayList<String> foldersToCombine,
													String sampleName) throws FileNotFoundException, IOException
	{
		ArrayList<BufferedReader> sampleReaders = new ArrayList<BufferedReader>();
		for (int f = 0; f < foldersToCombine.size(); f++)
		{
			BufferedReader fileReader = FileUtils.createReader(foldersToCombine.get(f) + sampleName);
			//fileReader.readLine();//skip header line
			if (new File(foldersToCombine.get(f) + sampleName).exists())
				sampleReaders.add(fileReader);
		}
		return sampleReaders;
	}

	private void executeSteps(	String matrixFn,
								String[] ourStudyStrings,
								String[] descriptions,
								ArrayList<String> folders,
								HashMap<String, HashMap<String, Integer>> sampleNameToFamilySampleNames,
								HashMap<String, String> patientSampleNameToFamilyId,
								HashMap<String, String> parentSampleNameToFamilyId,
								String rootFolder,
								String parentOutlierFolder,
								String childOutlierFolder,
								String maxReadCountOtherCasesFolder,
								String maxReadCountOtherParentFolder,
								HashMap<String, Double> maxBiosRawCount,
								HashMap<String, String> biosMedian) throws FileNotFoundException, IOException
	{
		if (matrixFn == null)
		{
			p("Warning matrixFn = null, skipping");
			return;
		}

		//get indexes of the different files
		StringFilters stringFiltersBios = new StringFilters();
		stringFiltersBios.addFilter(stringFiltersBios.new ForbiddenStringChecker(ourStudyStrings));
		int[] biosIndexes = getIndexes(	matrixFn,
										stringFiltersBios);
		StringFilters stringFiltersStudy = new StringFilters();
		stringFiltersStudy.addFilter(stringFiltersStudy.new RequiredStringChecker(ourStudyStrings));
		int[] studyIndexes = getIndexes(matrixFn,
										stringFiltersStudy);

		//get z-scores based Study samples from corrected data
		HashMap<String, int[]> sampleToOtherIndexesStudy = getStudyIndexesPerSample(matrixFn,
																					stringFiltersStudy,
																					sampleNameToFamilySampleNames);

		//get z-scores based BIOS samples from corrected data
		executeSubSteps(matrixFn,
						biosIndexes,
						studyIndexes,
						descriptions[0],
						folders.get(0),
						sampleToOtherIndexesStudy,
						false,
						patientSampleNameToFamilyId,
						parentSampleNameToFamilyId,
						rootFolder,
						sampleNameToFamilySampleNames,
						parentOutlierFolder,
						childOutlierFolder,
						maxReadCountOtherCasesFolder,
						maxReadCountOtherParentFolder,
						maxBiosRawCount,
						biosMedian);

		executeSubSteps(matrixFn,
						studyIndexes,
						studyIndexes,
						descriptions[1],
						folders.get(1),
						sampleToOtherIndexesStudy,
						true,
						null,
						null,
						rootFolder,
						sampleNameToFamilySampleNames,
						parentOutlierFolder,
						childOutlierFolder,
						maxReadCountOtherCasesFolder,
						maxReadCountOtherParentFolder,
						maxBiosRawCount,
						biosMedian);
	}

	private void executeSubSteps(	String matrixFn,
									int[] indexes,
									int[] writeIndexes,
									String description,
									String folder,
									HashMap<String, int[]> sampleToOtherIndexesStudy,
									boolean excludeFromReference,
									HashMap<String, String> patientSampleNameToFamilyId,
									HashMap<String, String> parentSampleNameToFamilyId,
									String rootFolder,
									HashMap<String, HashMap<String, Integer>> sampleNameToFamilySampleNames,
									String parentOutlierFolder,
									String casesOutlierFolder,
									String maxReadCountOtherCasesFolder,
									String maxReadCountOtherParentFolder,
									HashMap<String, Double> maxBiosRawCount,
									HashMap<String, String> biosMedian) throws FileNotFoundException, IOException
	{
		BufferedReader sampleReader = FileUtils.createReader(matrixFn);
		String writeFn = FileUtils.addBeforeExtention(	matrixFn,
														"_" + description + "Based_zscores");
		BufferedWriter zscoreWriter = FileUtils.createWriter(writeFn);
		String statsWriteFn = FileUtils.addBeforeExtention(	matrixFn,
															"_" + description + "Based_stats");
		BufferedWriter statsWriter = FileUtils.createWriter(statsWriteFn);

		String header = sampleReader.readLine();
		String[] colNames = header.split(	"\t",
											2)[1].split("\t");
		BufferedWriter[] studyWriters = createWriters(	writeIndexes,
														header,
														matrixFn,
														description,
														folder);

		writeSeparateHeaders(	writeIndexes,
								studyWriters,
								header,
								writeFn,
								description);

		zscoreWriter.write(header + "\n");
		statsWriter.write("RowName\tAverage\tStandard deviation\n");

		String line = null;
		HashMap<String, BufferedWriter> childOutlierWriters = null;
		HashMap<String, BufferedWriter> parentOutlierWriters = null;

		HashMap<String, BufferedWriter> maxReadCountOtherCasesWriters = null;
		HashMap<String, BufferedWriter> maxReadCountOtherParentWriters = null;

		while ((line = sampleReader.readLine()) != null)
		{
			HashMap<String, Integer> outlierCountOtherCases = createCountMatrixPatients(patientSampleNameToFamilyId);
			HashMap<String, Integer> outlierCountOtherParents = createCountMatrixPatients(patientSampleNameToFamilyId);

			HashMap<String, Integer> maxReadCountOtherCases = createCountMatrixPatients(patientSampleNameToFamilyId);
			HashMap<String, Integer> maxReadCountOtherParent = createCountMatrixPatients(patientSampleNameToFamilyId);
			//read the row
			Row row = Row.readRow(line);
			//get stdevs and avgs usign bios samples
			double[] average_Stdev_Largest_Median = row.getStatistics(	indexes,
																		false);
			double average = average_Stdev_Largest_Median[0];
			double stdev = average_Stdev_Largest_Median[1];
			double largest = average_Stdev_Largest_Median[2];
			double median = average_Stdev_Largest_Median[3];
			if (maxBiosRawCount != null)
				maxBiosRawCount.put(row.rowName,
									largest);

			if (biosMedian != null)
				biosMedian.put(	row.rowName,
								Double.toString(largest));

			Row rowZscores = new Row(	row.rowName,
										row.values.length);
			//for each sample 
			for (int c = 0; c < colNames.length; c++)
			{
				String sampleName = colNames[c];
				//exclude samples from the same family if calculating based on our study samples
				if (sampleToOtherIndexesStudy != null && excludeFromReference)
				{
					double[] average_Stdev = getNonFamilyStats(	row,
																sampleToOtherIndexesStudy,
																colNames,
																c);
					if (average_Stdev != null)
					{
						average = average_Stdev[0];
						stdev = average_Stdev[1];
					}
				}

				//z-score for this sample
				double zscore = (row.values[c] - average) / stdev;

				if (outlierCountOtherCases != null && patientSampleNameToFamilyId != null && patientSampleNameToFamilyId.containsKey(sampleName))
				{
					outlierCountOtherCases = countZscore(	outlierCountOtherCases,
															row.rowName,
															sampleName,
															zscore,
															patientSampleNameToFamilyId,
															sampleToOtherIndexesStudy,
															c,
															sampleNameToFamilySampleNames,
															this.zScoreCutoff);

					maxReadCountOtherCases = compareAndSetReadMaxCount(	maxReadCountOtherCases,
																		row.rowName,
																		sampleName,
																		(int) row.values[c],
																		patientSampleNameToFamilyId,
																		sampleToOtherIndexesStudy,
																		c,
																		sampleNameToFamilySampleNames,
																		this.zScoreCutoff);

				}

				if (outlierCountOtherParents != null && parentSampleNameToFamilyId != null && parentSampleNameToFamilyId.containsKey(sampleName))
				{
					outlierCountOtherParents = countZscore(	outlierCountOtherParents,
															row.rowName,
															sampleName,
															zscore,
															parentSampleNameToFamilyId,
															sampleToOtherIndexesStudy,
															c,
															sampleNameToFamilySampleNames,
															this.zScoreCutoff);

					maxReadCountOtherParent = compareAndSetReadMaxCount(outlierCountOtherParents,
																		row.rowName,
																		sampleName,
																		(int) row.values[c],
																		parentSampleNameToFamilyId,
																		sampleToOtherIndexesStudy,
																		c,
																		sampleNameToFamilySampleNames,
																		this.zScoreCutoff);

				}
				rowZscores.values[c] = zscore;
			}
			rowZscores.write(zscoreWriter);

			//write to separate files
			writeToSeparateFiles(	writeIndexes,
									studyWriters,
									rowZscores);
			//
			if (patientSampleNameToFamilyId != null)
			{
				new File(casesOutlierFolder).mkdirs();
				new File(parentOutlierFolder).mkdirs();

				if (maxReadCountOtherCasesFolder != null)
					new File(maxReadCountOtherCasesFolder).mkdirs();
				if (maxReadCountOtherParentFolder != null)
					new File(maxReadCountOtherParentFolder).mkdirs();

				childOutlierWriters = writeOutlierCounts(	outlierCountOtherCases,
															casesOutlierFolder,
															childOutlierWriters,
															rowZscores.rowName,
															matrixFn,
															"otherCases",
															description);
				parentOutlierWriters = writeOutlierCounts(	outlierCountOtherParents,
															parentOutlierFolder,
															parentOutlierWriters,
															rowZscores.rowName,
															matrixFn,
															"otherParents",
															description);

				if (maxReadCountOtherCasesFolder != null)
					maxReadCountOtherCasesWriters = writeOutlierCounts(	maxReadCountOtherCases,
																		maxReadCountOtherCasesFolder,
																		maxReadCountOtherCasesWriters,
																		rowZscores.rowName,
																		matrixFn,
																		"otherCases",
																		description);

				if (maxReadCountOtherParentFolder != null)
					maxReadCountOtherParentWriters = writeOutlierCounts(maxReadCountOtherParent,
																		maxReadCountOtherParentFolder,
																		maxReadCountOtherParentWriters,
																		rowZscores.rowName,
																		matrixFn,
																		"otherParents",
																		description);
			}
		}

		closeStreams(childOutlierWriters);
		closeStreams(parentOutlierWriters);
		closeStreams(studyWriters);
		statsWriter.close();
		zscoreWriter.close();
	}

	private HashMap<String, Integer> createCountMatrixPatients(HashMap<String, String> patientSampleNameToFamilyId)
	{
		if (patientSampleNameToFamilyId == null)
			return new HashMap<String, Integer>();

		HashMap<String, Integer> outlierCountOtherCases = new HashMap<String, Integer>();
		//add every patient to the list
		for (String patient : patientSampleNameToFamilyId.keySet())
			outlierCountOtherCases.put(	patient,
										0);

		return outlierCountOtherCases;
	}

	private void closeStreams(HashMap<String, BufferedWriter> writers) throws IOException
	{
		if (writers != null)
			for (Entry<String, BufferedWriter> writer : writers.entrySet())
			{
				writer.getValue().close();
			}
	}

	private HashMap<String, BufferedWriter> writeOutlierCounts(	HashMap<String, Integer> outlierCountOthers,
																String outlierFolder,
																HashMap<String, BufferedWriter> outlierWriters,
																String geneName,
																String matrixFn,
																String childOrParentString,
																String description) throws FileNotFoundException, IOException
	{
		if (outlierWriters == null)
		{
			outlierWriters = new HashMap<String, BufferedWriter>();
			int c = 0;
			String matrixName = FileUtils.removeExtention(new File(matrixFn).getName());
			for (String outlierCountOtherCase : outlierCountOthers.keySet())
			{
				BufferedWriter outlierCountOthersWriter = FileUtils.createWriter(outlierFolder + "/" + outlierCountOtherCase + ".txt");
				outlierWriters.put(	outlierCountOtherCase,
									outlierCountOthersWriter);
				outlierCountOthersWriter.write("geneOrTranscript-id\t" + matrixName + "_AbsolutezScoreAbove" + zScoreCutoff + "Counts_" + childOrParentString + "_" + description + "\n");
				c++;
			}
		}
		for (String outlierCountOther : outlierCountOthers.keySet())
		{
			int outliercount = outlierCountOthers.get(outlierCountOther);
			BufferedWriter outlierCountOtherCaseWriter = outlierWriters.get(outlierCountOther);
			outlierCountOtherCaseWriter.write(geneName + "\t" + outliercount + "\n");
		}

		return outlierWriters;

	}

	private HashMap<String, Integer> countZscore(	HashMap<String, Integer> outlierCountOthers,
													String gene,
													String sampleName,
													double zscore,
													HashMap<String, String> includeKeys,
													HashMap<String, int[]> sampleToOtherIndexesStudy,
													int c,
													HashMap<String, HashMap<String, Integer>> sampleNameToFamilySampleNames,
													double cutOff)
	{
		for (String sample : outlierCountOthers.keySet())
		{
			Integer count = outlierCountOthers.get(sample);
			if (count == null)
				count = new Integer(0);

			boolean isFamily = (sampleNameToFamilySampleNames.get(sample) != null && sampleNameToFamilySampleNames.get(sample).containsKey(sampleName));

			//if it is a studyReference sample that is not a family member and the outlier score is larger than 3
			if (Math.abs(zscore) >= cutOff && !isFamily && includeKeys.containsKey(sampleName))
			{
				count++;
			}
			outlierCountOthers.put(	sample,
									count);
		}
		return outlierCountOthers;
	}

	private HashMap<String, Integer> compareAndSetReadMaxCount(	HashMap<String, Integer> outlierCountOthers,
																String gene,
																String sampleName,
																int readCount,
																HashMap<String, String> includeKeys,
																HashMap<String, int[]> sampleToOtherIndexesStudy,
																int c,
																HashMap<String, HashMap<String, Integer>> sampleNameToFamilySampleNames,
																double cutOff)
	{
		for (String sample : outlierCountOthers.keySet())
		{
			Integer maxCount = outlierCountOthers.get(sample);
			if (maxCount == null)
				maxCount = new Integer(0);

			boolean isFamily = (sampleNameToFamilySampleNames.get(sample) != null && sampleNameToFamilySampleNames.get(sample).containsKey(sampleName));

			//if it is a studyReference sample that is not a family member and the outlier score is larger than 3
			if (Math.abs(readCount) > maxCount && !isFamily && includeKeys.containsKey(sampleName))
			{
				maxCount = readCount;
			}
			outlierCountOthers.put(	sample,
									maxCount);
		}
		return outlierCountOthers;
	}

	private double[] getNonFamilyStats(	Row row,
										HashMap<String, int[]> sampleNameToFamilySampleNames,
										String[] colNames,
										int c)
	{
		double[] average_Stdev = null;
		if (sampleNameToFamilySampleNames != null)
		{
			int[] nonFamilyIndexes = null;
			nonFamilyIndexes = sampleNameToFamilySampleNames.get(colNames[c]);

			if (nonFamilyIndexes != null)
			{
				average_Stdev = row.getStatistics(	nonFamilyIndexes,
													false);
			}
		}
		return average_Stdev;
	}

	private void writeSeparateHeaders(	int[] writeIndexes,
										BufferedWriter[] studyWriters,
										String header,
										String writeFn,
										String description) throws IOException
	{
		String colHeader = FileUtils.removeExtention(new File(writeFn).getName());

		String[] headers = header.split("\t");
		for (int i = 0; i < writeIndexes.length; i++)
		{
			int e = writeIndexes[i] + 1;
			studyWriters[i].write(headers[e].concat("\t").concat(colHeader).concat("\n"));
		}

	}

	private void writeToSeparateFiles(	int[] writeIndexes,
										BufferedWriter[] studyWriters,
										Row zScores) throws IOException
	{
		for (int i = 0; i < writeIndexes.length; i++)
		{
			int e = writeIndexes[i];
			studyWriters[i].write(zScores.getRowName().concat("\t").concat(Double.toString(zScores.getValues()[e])).concat("\n"));
		}

	}

	private void writeToSeparateFiles(	int[] writeIndexes,
										BufferedWriter[] studyWriters,
										Row zScores,
										String[] colNames,
										HashMap<String, ArrayList<String>> sampleNameToFamilyMemberSampleNames) throws IOException
	{
		for (int i = 0; i < writeIndexes.length; i++)
		{
			int e = writeIndexes[i];
			String sampleName = colNames[e];
			ArrayList<String> familyMemberSampleNames = sampleNameToFamilyMemberSampleNames.get(sampleName);

			studyWriters[i].write(zScores.getRowName().concat("\t").concat(Double.toString(zScores.getValues()[e])).concat("\n"));
		}

	}

	private BufferedWriter[] createWriters(	int[] writeIndexes,
											String header,
											String matrixFn,
											String description,
											String folder) throws FileNotFoundException, IOException
	{
		new File(folder).mkdir();

		BufferedWriter[] writers = new BufferedWriter[writeIndexes.length];
		String[] colNames = header.split(	"\t",
											2)[1].split("\t");
		for (int i = 0; i < writeIndexes.length; i++)
		{
			int c = writeIndexes[i];
			writers[i] = FileUtils.createWriter(folder + colNames[c] + ".txt");
		}

		return writers;
	}

	private void closeStreams(Closeable[] writers) throws IOException
	{
		for (Closeable writer : writers)
			writer.close();
	}

	private void closeStreams(ArrayList<BufferedReader> writers) throws IOException
	{
		for (Closeable writer : writers)
			writer.close();
	}

	private HashMap<String, int[]> getStudyIndexesPerSample(String correctedMatrixFn2,
															StringFilters stringFilters,
															HashMap<String, HashMap<String, Integer>> sampleNameToFamilySampleNames) throws FileNotFoundException, IOException
	{

		BufferedReader correctedReader = FileUtils.createReader(correctedMatrixFn);
		String header = correctedReader.readLine();
		String[] colNames = header.split(	"\t",
											2)[1].split("\t");//first split the first column off
		correctedReader.close();
		HashMap<String, int[]> studyIndexesPerSample = new HashMap<String, int[]>();

		for (int c = 0; c < colNames.length; c++)
		{
			if (stringFilters.check(colNames[c]))
			{
				String sampleName = colNames[c];
				HashMap<String, Integer> familySampleNames = sampleNameToFamilySampleNames.get(sampleName);

				if (familySampleNames != null)
				{
					int[] sampleIncludeIndexes = getIndexes(correctedMatrixFn2,
															stringFilters,
															familySampleNames);
					studyIndexesPerSample.put(	sampleName,
												sampleIncludeIndexes);
				}
				else//is probably not a child, but a parent
					;//p("Warning, sampleName: " + sampleName + " not found in sampleNameToFamilySampleNames file:\t" + this.sampleNameToFamilySampleNamesFn);

			}
		}

		return studyIndexesPerSample;
	}

	private int[] getIndexes(	String correctedMatrixFn2,
								StringFilters stringFilters) throws FileNotFoundException, IOException
	{
		return getIndexes(	correctedMatrixFn2,
							stringFilters,
							null);
	}

	private int[] getIndexes(	String correctedMatrixFn2,
								StringFilters stringFilters,
								HashMap<String, Integer> excludeSamples) throws FileNotFoundException, IOException
	{
		BufferedReader correctedReader = FileUtils.createReader(correctedMatrixFn);
		String header = correctedReader.readLine();
		String[] colNames = header.split(	"\t",
											2)[1].split("\t");//first split the first column off
		boolean[] include = new boolean[colNames.length];
		int n = 0;
		for (int c = 0; c < colNames.length; c++)
		{
			include[c] = stringFilters.check(colNames[c]);
			//exclude family samples			
			if (excludeSamples != null && excludeSamples.containsKey(colNames[c]))
				include[c] = false;

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

		correctedReader.close();
		return biosIndexes;
	}

	private List<Runnable> init()
	{
		//		//merge alle files samen die je wil gebruiken
		//		MergeFiles mergeFiles = new MergeFiles();
		//		mergeFiles.setFn1(correctedMatrixFn);
		//		mergeFiles.setFn1(fn2);
		//		mergeFiles.setKeepAllFromFn1(true);
		//		mergeFiles.setKeepAllFromFn1(false);

		List<Runnable> steps = new ArrayList<>();

		//		steps.add((Runnable) mergeFiles);

		return steps;
	}

	public double getMedian(double[] values)
	{
		Median median = new Median();
		double medianValue = median.evaluate(values);
		return medianValue;
	}
}
