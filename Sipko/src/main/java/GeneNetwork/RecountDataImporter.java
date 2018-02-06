package GeneNetwork;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import Tools.FileUtils;
import Tools.Script;

public class RecountDataImporter extends Script<ReactomeTermRetriever>
{
	//be sure to add the gtex and TGCA datasets as well as they are not on the sheet by default 
	String studyTablePath = "E:/Groningen/Data/GeneNetwork/Recount2/recount_selection_2018-01-16 12_35_54.csv";//from https://jhubiostatistics.shinyapps.io/recount/
	String baseDownloadUrlCountsV2 = "http://duffel.rail.bio/recount/v2/";
	String baseDownloadUrlPhenotypeFile = "http://duffel.rail.bio/recount/";
	
	String countsFileName = "counts_gene.tsv.gz";
	String infoFileName = "files_info.tsv";
	String writeFolderPerStudy = "E:/Groningen/Data/GeneNetwork/Recount2/studyDownloads/";
	
	@Override
	public void run()
	{
		try
		{
			if(new File(this.writeFolderPerStudy).exists())
				FileUtils.makeDir(writeFolderPerStudy);
			
			HashSet<String> studies = getStudyNamesFromStudyTable(this.studyTablePath);
//			
//			removeStudiesForWhichFolderAlreadyExists(studies);
//			//download counts files
//			downloadPerStudy(studies,this.countsFileName, true);
//			
//			studies = getStudyNamesFromStudyTable(this.studyTablePath);
//			//download info files for md5 checksums
//			downloadPerStudy(studies, this.infoFileName);
//			
//			//checkMd5checksums
//			removeStudiesWithMatchingMd5SumsFromHash(studies);
//			//////////////////////////////////////////////////////////////////////////////////////////////////////////
//			log("Second round of downloads; redownloading files that filed md5 checksum");
//			//second round of download counts files
//			downloadPerStudy(studies,this.countsFileName, true);
//			
//			studies = getStudyNamesFromStudyTable(this.studyTablePath);
//			//second round of download info files for md5 checksums
//			downloadPerStudy(studies, this.infoFileName);
//			//second round of checkMd5checksums
//			removeStudiesWithMatchingMd5SumsFromHash(studies);
			log("Downloading phenotype files");
			HashSet<String> studieMd5PhenotypeMd5_FailChecksum = downloadPerStudyPhenotypeFiles(studies);
			
//			if(studies.size()==0 && studieMd5PhenotypeMd5_FailChecksum.size()==0)
//				log("Done! Md5 checksums are matching for all files");
//			else
//			{
//				log("Done! Md5 checksums not matching for the following files (could also be files/URL are missing/bad):");
//				for(String study : studies)
//					log(study+"\t"+baseDownloadUrlCountsV2+study);
				log("Done! Md5 checksums not matching for the following phenotype/metric files (could also be files/URL are missing/bad):");
				for(String url: studieMd5PhenotypeMd5_FailChecksum)
				{
					log(url);
				}
//			}
			
			
			
		}catch(Exception e){e.printStackTrace();}
		log("Runtime of this script:\t" + this.getRunTime());
	}

	private HashSet<String> downloadPerStudyPhenotypeFiles(HashSet<String> studyNames) throws IOException
	{
		HashSet<String> studieMd5PhenotypeMd5_FailChecksum = new HashSet<String>();
		for(String studyName:studyNames)
		{
			String fileName=studyName+".tsv";
			String studyMatrixUrl =FileUtils.makeFolderNameEndWithSlash(baseDownloadUrlPhenotypeFile)+studyName+"/"+fileName;
			System.out.println("studyMatrixUrl=\t" + studyMatrixUrl);
			
			String writeFn = downloadFile(studyMatrixUrl, studyName, fileName);
			String infoFilePath = FileUtils.makeFolderNameEndWithSlash(writeFolderPerStudy)+infoFileName;
			
			boolean isPassChecksum=checkMd5CheckSum(writeFn,infoFilePath, fileName);
			if(!isPassChecksum)
			{
				log("Redownloading files. Md5 checksum was not matching for file:\t" + fileName);
				writeFn = downloadFile(studyMatrixUrl, studyName, fileName);
				String infoUrl=FileUtils.makeFolderNameEndWithSlash(baseDownloadUrlCountsV2)+studyName+"/"+infoFileName;
				log("infoUrl=\t" + infoUrl);
				infoFilePath = downloadFile(infoUrl, studyName, infoFileName);
				isPassChecksum=checkMd5CheckSum(writeFn,infoFilePath, fileName);
				
				//if they still dont match just print them at the end
				if(!isPassChecksum)
					studieMd5PhenotypeMd5_FailChecksum.add(infoUrl);
			}
		}
		return studieMd5PhenotypeMd5_FailChecksum;
	}

	private void removeStudiesWithMatchingMd5SumsFromHash(HashSet<String> studieNames) throws IOException
	{
		File[] folderPaths = new File(writeFolderPerStudy).listFiles();
		
		for(File folderPath: folderPaths)
		{
			if(!folderPath.isDirectory())
				continue;
			
			String countFilePath = folderPath+"/"+countsFileName;
			String infoFilePath = folderPath+"/"+infoFileName;
			
			boolean isPass=checkMd5CheckSum(countFilePath, infoFilePath, countsFileName);
			
			if(isPass)
			{
				String studyName = folderPath.getName();
				studieNames.remove(studyName);
			}
		}
		System.out.println("study_To_Urls.size= " + studieNames.size());
	}

	private boolean checkMd5CheckSum(	String countFilePath,
										String infoFilePath, String fileName) throws IOException
	{
		log("Checking md5 checksum for countFilePath"+ "\t" + infoFilePath);
		if(!new File(countFilePath).exists() || !new File(infoFilePath).exists() || new File(countFilePath).length()==0 || new File(infoFilePath).length()==0)
			return false;
		
		String countFileChecksum = FileUtils.getMd5(countFilePath);
		String l = FileUtils.getLine(infoFilePath, fileName);
		String md5Checksum = l.split("\t")[1];
		
		BufferedWriter checkSumWriter = FileUtils.createWriter(countFilePath+".md5");
		checkSumWriter.write(countFileChecksum);
		checkSumWriter.close();
		
		if(md5Checksum.equals(countFileChecksum))
			return true;
		return false;
	}

	private void removeStudiesForWhichFolderAlreadyExists(HashSet<String> studyNames)
	{
		ArrayList<String> studiesToRemove = new ArrayList<String>();
		for(String studyName:studyNames)
		{
			if(new File(this.writeFolderPerStudy+studyName).exists())
			{
				studiesToRemove.add(studyName);
			}
		}
		for(String studyName: studiesToRemove)
			studyNames.remove(studyName);
		
	}

	private void downloadPerStudy(HashSet<String> study_To_Urls, String fileName) throws IOException
	{
		downloadPerStudy(study_To_Urls, fileName, false);
	}
	private void downloadPerStudy(HashSet<String> studyNames, String fileName, boolean verbose) throws IOException
	{
		FileUtils.makeDir(this.writeFolderPerStudy);
		int n =1;
		for(String studyName:studyNames)
		{
			if(verbose)
				log("Downloading study: "+ n +"\t" + studyName);
			String countMatrixUrl =FileUtils.makeFolderNameEndWithSlash(baseDownloadUrlCountsV2)+studyName+"/"+fileName;
			
			downloadFile(countMatrixUrl, studyName, fileName);
			
			n++;
		}		
	}

	private String downloadFile(String url, String studyName, String fileName)
	{
		log(url);
		String writeFolder = this.writeFolderPerStudy+studyName+"/";//do not change this, script assumes studyname is the upper foldername to do the md5 checksums
		FileUtils.makeDir(writeFolder);
		String writeFn = writeFolder+fileName;
		try
		{
			FileUtils.saveFileFromUrl(url, writeFn);
		}catch(Exception e)
		{
			e.printStackTrace();
			log("Skipping study due to Invalid URL:\t" + url);
		}
		return writeFn;
	}

	private HashSet<String> getStudyNamesFromStudyTable(String studyTablePath) throws IOException
	{
		BufferedReader studyTableReader = FileUtils.createReader(studyTablePath);
		String header = studyTableReader.readLine();
		String line = null;
		HashSet<String> studies = new HashSet<String>();
		
		while((line=studyTableReader.readLine())!=null)
		{
			String studyName = line.split(",")[0].replace("\"", "");
			studies.add(studyName);
		}
		
		return studies;
	}
}

