package TextEditing;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;

import Tools.FileSearcher;
import Tools.FileUtils;
import Tools.Script;

public class FastQCmerger extends Script<FastQCmerger>
{
	String fastQcFolder = null;
	String writeFn = null;

	@Override
	public void run()
	{
		try
		{
			init();
			
			BufferedWriter duplicationSummary = FileUtils.createWriter(writeFn);
			duplicationSummary.write("SampleName\tTotal Deduplicated Percentage\n");
			FileSearcher fastQcSearcher = new FileSearcher();
			fastQcSearcher.setSearchStrings(new String[]{"fastqc_data.txt"});
			fastQcSearcher.setFolders(fastQcFolder);
			fastQcSearcher.run();
			
			
			String fastQcFileFns = fastQcSearcher.getWriteName();
			p(fastQcFileFns);			
			
			BufferedReader fastQcFilesReader = FileUtils.createReader(fastQcFileFns);
			String fastqFn = null;
			while((fastqFn = fastQcFilesReader.readLine())!=null)
			{
				BufferedReader fastQcFileReader = FileUtils.createReader(fastqFn);
				String line = null;
				while((line =fastQcFileReader.readLine())!= null)
				{
					if(line.startsWith("#Total Deduplicated Percentage"))
					{
						String duplicationPercentage = line.split("\t")[1];
						
						duplicationSummary.write(FileUtils.getFolderName(fastqFn)+"\t"+duplicationPercentage+"\n");
					}
				}
			}
			fastQcFilesReader.close();
			duplicationSummary.close();
		}catch(Exception e){e.printStackTrace();}
		
	}

	private void init()
	{
		if(writeFn ==null)
			writeFn = fastQcFolder+"duplicationSummary.txt";
		
	}
}
