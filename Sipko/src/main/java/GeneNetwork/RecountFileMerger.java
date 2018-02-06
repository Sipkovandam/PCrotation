package GeneNetwork;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.HashMap;
import java.util.HashSet;

import MatrixScripts.MyMatrix;
import Tools.FileUtils;
import Tools.Script;

public class RecountFileMerger extends Script<RecountFileMerger>
{
	String folderName = null;
	String countFileName = null;
	String writeFn = null;
	String writeFnMetricMatrix = null;
	String metrixHeaderToUseFn = null;
	
	
	@Override
	public void run()
	{
		try
		{
			File[] recountFolderPaths = new File(folderName).listFiles();
//			
//			MyMatrix totalMatrix = null;
//
//			int n =0;
//			
//			//merge count matrixes
//			for(File recountFolderPath : recountFolderPaths)
//			{
//				n++;
//				if(!recountFolderPath.isDirectory())
//					continue;
//				String fileName = FileUtils.makeFolderNameEndWithSlash(recountFolderPath.getAbsolutePath())+countFileName;
//				if(!new File(fileName).exists() || new File(fileName).length()==0)
//					continue;
//				
//				log("Merging file\t" + n + "/" + recountFolderPaths.length + "\t" + fileName);
//				if(totalMatrix==null)
//					totalMatrix=new MyMatrix(fileName, -1);
//				else
//				{
//					MyMatrix countMatrix = new MyMatrix(fileName, -1);
//					totalMatrix=totalMatrix.mergeColumns(countMatrix);
//				}
//			}
//			log("Writing totalMatrix:\t" + writeFn);
//			totalMatrix.write(writeFn);
			
			//merge phenotype/metric matrixes
			String headerLine = FileUtils.getLine(metrixHeaderToUseFn, "\t");
			String[] colNames = headerLine.split("\t");
			
			BufferedWriter totalMetricMatrix = FileUtils.createWriter(writeFnMetricMatrix);
			totalMetricMatrix.write(headerLine+"\n");
			HashSet<String> missingMetricFiles = new HashSet<String>();
			for(File recountFolderPath : recountFolderPaths)
			{
				
				if(!recountFolderPath.isDirectory())
					continue;
				String sampleName = recountFolderPath.getName();
				String fileName = FileUtils.makeFolderNameEndWithSlash(recountFolderPath.getAbsolutePath())+sampleName+".tsv";
				if(!new File(fileName).exists())
				{
					missingMetricFiles.add(fileName);
					continue;
				}
				
				log("File=\t" + fileName);
				BufferedReader reader = FileUtils.createReader(fileName);
				String line = reader.readLine();
				HashMap<String, Integer> colName_To_Colnubmer=FileUtils.makeHash(line);
				while((line=reader.readLine())!=null)
				{
					if(line.startsWith("\""))//in some files the format was corrupted. E.g. E:\Groningen\Data\GeneNetwork\Recount2\studyDownloads\SRP029334/
					{
						log("In this file:\t" + fileName +"\t this corrupt line was detected:\t" + line + "\tin file\t");
						continue;
					}
					String[] colValues=line.split("\t");
					int c =0;
					for(String colName: colNames)
					{
						if(!colName_To_Colnubmer.containsKey(colName))
							continue;
						int col = colName_To_Colnubmer.get(colName);
						if(c!=0)
							totalMetricMatrix.write("\t");
						
						if(col<colValues.length)
							totalMetricMatrix.write(colValues[col]);
						else
							log("In this file:\t" + fileName +"\t this corrupt line was detected:\t" + line + "\tin file\t");
						
						c++;
					}
					if(c!=0)
						totalMetricMatrix.write("\n");
				}
			}
			if(missingMetricFiles.size()>0)
			{
				log("The following metric files are missing:");
				for(String missingMetricFile : missingMetricFiles)
					log(missingMetricFile);
			}
			totalMetricMatrix.close();	
		}catch(Exception e)
		{
			e.printStackTrace();
		}
	}
}
