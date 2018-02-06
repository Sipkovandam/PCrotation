package PCA;

import java.io.BufferedReader;
import java.io.BufferedWriter;

import MatrixScripts.MyMatrix;
import Tools.FileUtils;
import Tools.Script;
import htsjdk.samtools.cram.encoding.writer.Writer;

public class RemoveNoVarianceRows extends Script<RemoveNoVarianceRows>
{
	//removes rows that have no variance. Assumes matrix has rownames and colnames.
	String fileName = null;
	String writeName = null;
	
	@Override
	public void run()
	{
		try
		{
			BufferedReader fileReader=FileUtils.createReader(fileName);
			BufferedWriter fileWriter=FileUtils.createWriter(writeName);
			
			String line = fileReader.readLine();
			fileWriter.write(line);
			log("Removing rows without variance");
			int nNoVariance = 0;
			while((line=fileReader.readLine())!=null)
			{
				String[] eles = line.split("\t");
				
				boolean hasVariance=isVariance(eles);
				
				if(hasVariance)
					fileWriter.write(line+"\n");
				else
					nNoVariance++;
			}
			fileWriter.close();
			log("Removed this many rows (these had no variance): \t" + nNoVariance);
		}catch(Exception e){e.printStackTrace();}
	}

	private boolean isVariance(String[] eles)
	{
		for(int c=2; c < eles.length;c++)
		{
			if(!eles[1].equals(eles[c]))
				return true;
		}
		return false;
	}
}
