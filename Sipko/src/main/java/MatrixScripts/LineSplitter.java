package MatrixScripts;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;

import Tools.FileUtils;
import Tools.Script;

public class LineSplitter extends Script<LineSplitter>
{
	//splits rows based on a particular column with ids separated by ",". 
	//Can also get the names to split by from a separate file having the splice names in the first column and the comma separated gene name list in the second (supply <conversionFn>)
	
	String fn = null;
	String conversionFn = null;
	int columnWithIds = 0;
	String separator = ",";
	String writeFn = null;
	
	@Override
	public void run()
	{
		try
		{
			HashMap<String,String> spliceToGenename = null;
			//hash that gets the genes that this splice juction overlaps
			if(conversionFn != null)
				spliceToGenename=FileUtils.readStringStringHash(conversionFn);
				
			BufferedWriter writer = FileUtils.createWriter(writeFn);
			BufferedReader reader = FileUtils.createReader(fn);
			String line = reader.readLine();
			writer.write("geneName\t"+line+"\n");
			while((line=reader.readLine())!=null)
			{
				String spliceName= line.split("\t")[0];
				String[] geneIds = null;
				if(conversionFn==null)
					geneIds=line.split("\t")[columnWithIds].split(separator);
				else
				{
					if(spliceToGenename.get(spliceName)!=null)
						geneIds=spliceToGenename.get(spliceName).replace("\"", "").split(separator);
					else
						geneIds=new String[]{"null"};
				}
				for(String geneId : geneIds)
					writer.write(geneId+"\t"+line+"\n");
			}
			
			writer.close();
			reader.close();
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
}
