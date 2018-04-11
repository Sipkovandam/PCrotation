package Gdio;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import Tools.ExcelFileParser;
import Tools.FileUtils;
import Tools.Script;

public class MolgenisEmxParser extends Script<MolgenisEmxParser>
{
	//Script has special treatment for ...EntryDate :)\
	//still need to remove duplicate values from persons tabel after running this script
	//assumes forwardread files end with: _1.fq.gz and backward read files end with _2.fq.gz
	String sipkoFormattedSheet = "";
	String emptyEmxFile = "";
	String writeFn = "";
	String packageName = "gdio_";

	transient HashMap<String, String> dnaNumber_To_PersonId = new HashMap<String, String>();
	transient HashMap<String, String> idName_To_SheetName = new HashMap<String, String>();
	transient HashMap<String, String> collection_To_CollectionId = new HashMap<String, String>();
	transient HashMap<String, HashMap<String, Integer>> table_To_ColNames_To_ColNumber = new HashMap<String, HashMap<String, Integer>>();
	transient HashMap<String, String> fileName_To_FileID = new HashMap<String, String>();

	transient HashMap<String, Integer> columnID_To_currentNumber = new HashMap<String, Integer>();
	transient HashSet<String> alreadyAddedLines = new HashSet<String>();

	@Override
	public void run()
	{
		try
		{
			log("Started parsing file");
			BufferedReader sipkoFormattedSheetReader = FileUtils.createReader(sipkoFormattedSheet);
			ExcelFileParser excelFileWriter = new ExcelFileParser();
			excelFileWriter.createWriter(writeFn);

			parseEmptyEmxFile(	this.emptyEmxFile,
								excelFileWriter);

			String header = sipkoFormattedSheetReader.readLine();
			String line = null;
			while ((line = sipkoFormattedSheetReader.readLine()) != null)
			{
				String[] eles = line.replace(	"\"",
												"").split("\t");
				if (eles.length > 2)
				{
					addPersonDnaNumberToPersonId(eles);
				}

			}

			log("Persons=" + dnaNumber_To_PersonId.size());
			sipkoFormattedSheetReader = FileUtils.createReader(sipkoFormattedSheet);
			header = sipkoFormattedSheetReader.readLine();
			while ((line = sipkoFormattedSheetReader.readLine()) != null)
			{
				String[] cells = line.replace(	"\"",
												"").split("\t");

				HashMap<String, String> featureName_To_Feature = new HashMap<String, String>();
				addFeature(	featureName_To_Feature,
							cells,
							"dnaNumber",
							0);
				addFeature(	featureName_To_Feature,
							cells,
							"OriginalID",
							1);
				addFeature(	featureName_To_Feature,
							cells,
							"FileColumnSampleId",
							1);
				addFeature(	featureName_To_Feature,
							cells,
							"fastqsForward",
							2);
				addFeature(	featureName_To_Feature,
							cells,
							"fastqsBackward",
							4);
				addFeature(	featureName_To_Feature,
							cells,
							"MaterialType",
							6);
				addFeature(	featureName_To_Feature,
							cells,
							"HpoCode",
							7);
				addFeature(	featureName_To_Feature,
							cells,
							"Mutation",
							8);
				addFeature(	featureName_To_Feature,
							cells,
							"father",
							9);
				addFeature(	featureName_To_Feature,
							cells,
							"mother",
							10);
				addFeature(	featureName_To_Feature,
							cells,
							"FamilyID",
							11);
				addFeature(	featureName_To_Feature,
							cells,
							"Sex",
							12);
				addFeature(	featureName_To_Feature,
							cells,
							"TissueType",
							13);
				addFeature(	featureName_To_Feature,
							cells,
							"CollectionName",
							14);
				addFeature(	featureName_To_Feature,
							cells,
							"connectedFiles1",
							17);
				addFeature(	featureName_To_Feature,
							cells,
							"connectedFiles1_SampleName",
							19);

				//skip lines for which no fastq files exist
				if (featureName_To_Feature.get("fastqsForward") == null || featureName_To_Feature.get("fastqsForward").equals(""))
					continue;

				//"dd-MM-yyyy"
				DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy-MM-dd");
				LocalDate localDate = LocalDate.now();
				featureName_To_Feature.put(	"EntryDate",
											dtf.format(localDate));

				featureName_To_Feature.put(	"PersonID",
											dnaNumber_To_PersonId.get(featureName_To_Feature.get("dnaNumber")));
				featureName_To_Feature.put(	"PseudoID",
											dnaNumber_To_PersonId.get(featureName_To_Feature.get("dnaNumber")));
				featureName_To_Feature.put(	"MotherID",
											dnaNumber_To_PersonId.get(featureName_To_Feature.get("mother")));
				
				featureName_To_Feature.put(	"FatherID",
											dnaNumber_To_PersonId.get(featureName_To_Feature.get("father")));

				String collectionId = collection_To_CollectionId.get(featureName_To_Feature.get("CollectionName"));
				if (collectionId == null)
				{
					collectionId = getNewIDnumber("CollectionID");
					collection_To_CollectionId.put(	featureName_To_Feature.get("CollectionName"),
													collectionId);
				}

				featureName_To_Feature.put(	"CollectionID",
											collectionId);
				featureName_To_Feature.put(	"ConsentID",
											"");

				//add person to personSheet
				addToSheet(	"gdio_Person",
							featureName_To_Feature,
							excelFileWriter);
				addToSheet(	"gdio_Sample",
							featureName_To_Feature,
							excelFileWriter);
				addToSheet(	"gdio_Collection",
							featureName_To_Feature,
							excelFileWriter);
				addToSheet(	"gdio_CollectionPerson",
							featureName_To_Feature,
							excelFileWriter);

				if (hasValue(featureName_To_Feature.get("Mutation")) || hasValue(featureName_To_Feature.get("HpoCode")))
				{
					featureName_To_Feature.put(	"AssessmentID",
					                           	getNewIDnumber("AssessmentID"));
					featureName_To_Feature.put(	"gdio_Diagnosis",
					                           	getNewIDnumber("gdio_Diagnosis"));

					addToSheet(	"gdio_Assessment",
								featureName_To_Feature,
								excelFileWriter);
					
					if (hasValue(featureName_To_Feature.get("Mutation")))
					{
						addToSheet(	"gdio_Mutation",
									featureName_To_Feature,
									excelFileWriter);

						addToSheet(	"gdio_Diagnosis",
									featureName_To_Feature,
									excelFileWriter);
					}

					if (featureName_To_Feature.get("HpoCode") != null)
					{
						String[] hpoCodes = featureName_To_Feature.get("HpoCode").split(",");
						addToSheet(	"gdio_Assessment",
										featureName_To_Feature,
										excelFileWriter);
						
						for (String hpoCode : hpoCodes)
						{
							featureName_To_Feature.put(	"HpoCode",
														hpoCode);
							featureName_To_Feature.put(	"HpoID",
							                           	getNewIDnumber("HpoID"));
							featureName_To_Feature.put(	"HpoAssessmentID",
							                           	getNewIDnumber("HpoAssessmentID"));
							
							addToSheet(	"gdio_Hpo",
										featureName_To_Feature,
										excelFileWriter);
							addToSheet(	"gdio_HpoAssessment",
										featureName_To_Feature,
										excelFileWriter);
						}			
					}
				}
				

				//add files to filesheet
				String[] forwardFiles = featureName_To_Feature.get("fastqsForward").split(",");
				
				boolean newSample = true; 
				for (String forwardFile : forwardFiles)
				{
					featureName_To_Feature.put(	"FileColumnSampleID",
												"");
					featureName_To_Feature.put(	"FileName",
												forwardFile);
					featureName_To_Feature.put(	"FileID",
												getNewIDnumber("FileID"));
					featureName_To_Feature.put(	"SampleFileID",
												getNewIDnumber("SampleFileID"));
					String fileTypeForwardFile = getFileType(forwardFile);
					featureName_To_Feature.put(	"FileType",
					                           	fileTypeForwardFile);
										
					addToSheet(	"gdio_File",
								featureName_To_Feature,
								excelFileWriter);
					addToSheet(	"gdio_SampleFile",
								featureName_To_Feature,
								excelFileWriter);

					addConnectedFiles(	featureName_To_Feature.get("FileID"),
										featureName_To_Feature.get("connectedFiles1"),
										featureName_To_Feature.get("connectedFiles1_SampleName"),
										featureName_To_Feature,
										excelFileWriter, newSample);
					//backwardReads
					if (forwardFile.endsWith("_1.fq.gz"))
					{
						String backwardName=forwardFile.replace("_1.fq.gz",
								"_2.fq.gz");
						featureName_To_Feature.put(	"FileColumnSampleID",
								"");
						featureName_To_Feature.put(	"FileName",backwardName
													);
						String fileTypeBackward = getFileType(backwardName);
						featureName_To_Feature.put(	"FileType",
						                           	fileTypeBackward);
						
						featureName_To_Feature.put(	"SampleFileID",
													getNewIDnumber("SampleFileID"));
						
						featureName_To_Feature.put(	"ParentFileID",
													featureName_To_Feature.get("FileID"));
						featureName_To_Feature.put(	"FileID",
													getNewIDnumber("FileID"));
						featureName_To_Feature.put(	"ChildFileID",
													featureName_To_Feature.get("FileID"));
						featureName_To_Feature.put(	"FileRelationID",
													getNewIDnumber("FileRelationID"));
						addToSheet(	"gdio_File",
									featureName_To_Feature,
									excelFileWriter);
						addToSheet(	"gdio_FileRelation",
									featureName_To_Feature,
									excelFileWriter);
						addToSheet(	"gdio_SampleFile",
									featureName_To_Feature,
									excelFileWriter);

						addConnectedFiles(	featureName_To_Feature.get("FileID"),
											featureName_To_Feature.get("connectedFiles1"),
											featureName_To_Feature.get("connectedFiles1_SampleName"),
											featureName_To_Feature,
											excelFileWriter, false);
					}
					newSample=false;
				}
				addToSheet(	"gdio_AnonymizerHash",
							featureName_To_Feature,
							excelFileWriter);

				addToSheet(	"gdio_AnonymizerHash",
							featureName_To_Feature,
							excelFileWriter);

			}
			excelFileWriter.closeWriter();
			sipkoFormattedSheetReader.close();
		} catch (FileNotFoundException e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	private String getFileType(String fileName)
	{
		String fileType= fileName.replace(".gz", "").replaceAll(".*\\.","");
		if(fileType.equals("fq"))
			fileType="fastq";
		return fileType;
	}

	private void addConnectedFiles(	String parentFileID,
									String childFn,
									String sampleNameInFile,
									HashMap<String, String> featureName_To_Feature,
									ExcelFileParser excelFileWriter, boolean newSample)
	{
		if(childFn==null || childFn.equals(""))
			return;
		
		HashMap<String, String> featureName_To_Feature2 = (HashMap<String, String>) featureName_To_Feature.clone();
		
		String fileType = getFileType(childFn);
		featureName_To_Feature2.put(	"FileType",
		                           	fileType);
		
		featureName_To_Feature2.put(	"FileRelationID",
									getNewIDnumber("FileRelationID"));
		featureName_To_Feature2.put(	"ParentFileID",
		                           	parentFileID);
		
		
		featureName_To_Feature2.put(	"FileColumnSampleID",
									sampleNameInFile);
		featureName_To_Feature2.put(	"FileName",
		                           	childFn);
		
		//important this one is done last // this also adds new files to the file table if needed
		String connectedFileID=getConnectedFileID(childFn,
		                                          featureName_To_Feature2,
			                    					excelFileWriter);
		
		featureName_To_Feature2.put(	"FileID",
			                           	connectedFileID);	
		
		featureName_To_Feature2.put(	"ChildFileID",
		                           	connectedFileID);
		
		if(newSample)
		{
			featureName_To_Feature2.put(	"SampleFileID",
										getNewIDnumber("SampleFileID"));
			
			addToSheet(	"gdio_SampleFile",
			           	featureName_To_Feature2,
						excelFileWriter);
		}
		
		addToSheet(	"gdio_FileRelation",
		           	featureName_To_Feature2,
					excelFileWriter);
	}

	private String getConnectedFileID(String childFn, HashMap<String, String> featureName_To_Feature, ExcelFileParser excelFileWriter)
	{
		String connectedFileID= fileName_To_FileID.get(childFn);
		if(connectedFileID==null)
		{
			connectedFileID=getNewIDnumber("FileID");
			featureName_To_Feature.put(	"FileID",
			                           	connectedFileID);
			
			addToSheet(	"gdio_File",
						featureName_To_Feature,
						excelFileWriter);
			

			fileName_To_FileID.put(childFn, connectedFileID);
		}
		return connectedFileID;
	}

	private boolean hasValue(String string)
	{
		if (string != null && !string.equals(""))
			return true;
		return false;
	}

	private void addPersonDnaNumberToPersonId(String[] cells)
	{
		String dnaNumberOrRnaNumber = cells[0];
		String forwardFastqs = cells[2];
		if (this.dnaNumber_To_PersonId.containsKey(dnaNumberOrRnaNumber) || forwardFastqs == null || forwardFastqs.equals(""))
			return;
		String newIdForColumn = getNewIDnumber("PersonID");
		this.dnaNumber_To_PersonId.put(	dnaNumberOrRnaNumber,
										newIdForColumn);
	}

	private String getId(	HashMap<String, String> dnaNumber_To_PersonId2,
							HashMap<String, String> featureName_To_Feature,
							ExcelFileParser excelFileWriter,
							String dnaNumber,
							String column,
							String[] sheets)
	{
		if (dnaNumber_To_PersonId2.containsKey(dnaNumber))
			return dnaNumber_To_PersonId2.get(dnaNumber);
		String newIdForColumn = getNewIDnumber(column);
		dnaNumber_To_PersonId2.put(	dnaNumber,
									newIdForColumn);

		return newIdForColumn;
	}

	private void addFeature(HashMap<String, String> featureName_To_Feature,
							String[] cells,
							String colName,
							int col)
	{
		if (col < cells.length && cells[col].length() > 0)
			featureName_To_Feature.put(	colName,
										cells[col]);
	}

	private void addToSheet(String gdioSheetName,
							HashMap<String, String> featureName_To_Feature,
							ExcelFileParser excelFileWriter)
	{
		HashMap<String, Integer> colNames_To_ColNumber = table_To_ColNames_To_ColNumber.get(gdioSheetName);

		String line = insertEachColumn(	colNames_To_ColNumber,
										featureName_To_Feature);

		if (line != null && !alreadyAddedLines.contains(line))
		{
			excelFileWriter.write(	line + "\n",
									gdioSheetName);
			alreadyAddedLines.add(line);
		}
	}

	private String insertEachColumn(HashMap<String, Integer> colNames_To_ColNumber,
									HashMap<String, String> featureName_To_Feature)
	{
		Set<String> colNames = colNames_To_ColNumber.keySet();
		String[] colValues = new String[colNames.size()];

		if (!isCanAddValuesTosheet(	featureName_To_Feature,
									colNames))
			return null;

		for (String colName : colNames)
		{
			if (!featureName_To_Feature.containsKey(colName))
				if (colName.contains("ID") && !colName.contains("ensemblID") && !colName.contains("FamilyID"))
				{
					String newColIDnumber = getNewIDnumber(colName);
					featureName_To_Feature.put(	colName,
												newColIDnumber);
				}
				else
					continue;

			int colNumber = colNames_To_ColNumber.get(colName);
			colValues[colNumber] = featureName_To_Feature.get(colName);
		}

		String line = FileUtils.makeLineFromArray(colValues);
		return line;
	}

	private boolean isCanAddValuesTosheet(	HashMap<String, String> featureName_To_Feature,
											Set<String> colNames)
	{
		for (String colName : colNames)
		{
			if (featureName_To_Feature.containsKey(colName) && !colName.contains("EntryDate"))
			{
				return true;
			}
		}
		return false;
	}

	private String getNewIDnumber(String colName)
	{
		int currentNumber = 0;
		if (columnID_To_currentNumber.containsKey(colName))
			currentNumber = columnID_To_currentNumber.get(colName);

		columnID_To_currentNumber.put(	colName,
										currentNumber + 1);
		return colName + "_" + currentNumber;
	}

	private void parseEmptyEmxFile(	String emptyEmxFile,
									ExcelFileParser excelFileWriter) throws IOException
	{
		ExcelFileParser emptyEmxReader = new ExcelFileParser();
		ArrayList<String> sheetNames = emptyEmxReader.getSheetNames(emptyEmxFile);

		for (int s = 0; s < sheetNames.size(); s++)
		{
			String sheetName = sheetNames.get(s);
			if (sheetName.contains(packageName))
				continue;

			parseSheet(	emptyEmxFile,
						sheetName,
						s,
						excelFileWriter);
		}
	}

	private void parseSheet(String emptyEmxFile2,
							String sheetName,
							int sheetNumber,
							ExcelFileParser excelFileWriter) throws IOException
	{
		ExcelFileParser emptyEmxReader = new ExcelFileParser();
		emptyEmxReader.getReader(	emptyEmxFile2,
									sheetNumber);
		String line = emptyEmxReader.readLine();
		String[] columns = line.split("\t");
		excelFileWriter.write(	line + "\n",
								sheetName);

		HashMap<String, Integer> ColNames_To_ColNumber = new HashMap<String, Integer>();
		for (int c = 0; c < columns.length; c++)
		{
			ColNames_To_ColNumber.put(	columns[c],
										c);
			if (columns[c].contains("EntryDate"))//special treatment for entry date :)
				ColNames_To_ColNumber.put(	"EntryDate",
											c);
		}
		//fully copy first 3 sheets
		if (!sheetName.contains(packageName))
		{
			while ((line = emptyEmxReader.readLine()) != null && line.replace(	"\t",
																				"").length() > 0)
			{
				excelFileWriter.write(	line + "\n",
										sheetName);
			}
		}

		table_To_ColNames_To_ColNumber.put(	sheetName,
											ColNames_To_ColNumber);
		log("SheetName=\t" + sheetName + "\tsheetNumber=\t" + sheetNumber + "\tline =\t" + line);
	}
}
