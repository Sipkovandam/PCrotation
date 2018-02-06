package Kallisto;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

import Tools.FileSearcher;
import Tools.FileUtils;
import Tools.Script;

public class MappingPercentagesRetriever extends Script<MappingPercentagesRetriever>
{
	String kallistoFolder = null;
	String writeFn = null;
	
	@Override
	public void run()
	{
		try
		{
			if(writeFn==null)
				writeFn=FileUtils.makeFolderNameEndWithSlash(kallistoFolder) + "mappingPercentages.txt";
			
			String errorFilesFn = FileUtils.makeFolderNameEndWithSlash(kallistoFolder) + "errorFileNames.txt";
			FileSearcher fileSearcher = new FileSearcher();
			fileSearcher.setFolders(kallistoFolder);
			fileSearcher.setSearchStrings(new String[]{".err"});
			fileSearcher.setWriteName(errorFilesFn);
			
			fileSearcher.run();
			
			BufferedReader fnReader = FileUtils.createReader(errorFilesFn);
			String line = null;
			BufferedWriter percentageWriter = FileUtils.createWriter(writeFn);
			percentageWriter.write("FileName\tMappedReads\tTotalReads\tPercentage\n");
			while((line = fnReader.readLine())!= null)
			{
				if(!new File(line).exists())
					continue;
				

				BufferedReader errorReader = FileUtils.createReader(line);
				
				//get to the line that has the mapping percentages
				String errorLine = errorReader.readLine();
				while(errorLine!=null && !errorLine.contains(" pseudoaligned"))
					errorLine=errorReader.readLine();
				
				if(errorLine==null)
				{
					log("Warning, this file does not contain mapping percentages:\n" + line);
					continue;
				}

				String total = null;
				String mapped = null;
				double percentage = 0;
				if(errorLine.contains("[quant] processed "))
				{
					total = errorLine.split("processed ")[1].split(" reads, ")[0].replace(",", "");
					mapped = errorLine.split(" reads, ")[1].split(" reads ")[0].replace(",", "");
					percentage = Double.parseDouble(mapped)/Double.parseDouble(total);
				}
				percentageWriter.write(line+"\t" + mapped+"\t" + total + "\t"+ percentage+"\n");
			}
			
			percentageWriter.close();
		}catch(Exception e){e.printStackTrace();}
		
	}
}
