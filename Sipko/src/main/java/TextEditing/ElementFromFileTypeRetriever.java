package TextEditing;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;

import Tools.FileSearcher;
import Tools.FileUtils;
import Tools.Script;

public class ElementFromFileTypeRetriever extends Script<ElementFromFileTypeRetriever>
{
	String lineMustContain = null;
	int column = -1;
	String writeFn = null;
	String header = "SampleName\tMetric";
	
	FileSearcher fileSearcher = new FileSearcher();
	
	@Override
	public void run()
	{
		try
		{
			init();
		
			fileSearcher.run();

			if(writeFn == null)
				writeFn=FileUtils.makeFolderNameEndWithSlash(fileSearcher.getFolders())+"summary.txt";
			
			ArrayList<String> fns = FileUtils.readArrayList(fileSearcher.getWriteName());
			BufferedWriter writer = FileUtils.createWriter(writeFn);
			writer.write(this.header+"\n");
			for(String fn : fns)
			{
				BufferedReader reader = FileUtils.createReader(fn);
				String line = null;
				String retrievedElement = null;
				while((line = reader.readLine()) != null)
				{
					
					if(lineMustContain!=null && !line.contains(lineMustContain))
						continue;
					System.out.println(line);
					if(column >0)
						retrievedElement=line.split("\t")[column];
					
					System.out.println(retrievedElement);
				}
				writer.write(FileUtils.getFolderName(fn)+"\t"+retrievedElement+"\n");
				reader.close();
			}
			writer.close();
		}catch(Exception e){e.printStackTrace();}
	}

	private void init()
	{
	}
}
