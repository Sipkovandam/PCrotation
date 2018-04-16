package Tests;

import java.io.BufferedReader;

import Tools.FileUtils;

public class PrintLargeFile
{
	public static void main(String[] args)
	{
		try
		{
			String fn = "E:/Groningen/Data/Annotation/GRCh38/GeneNetwork/hpo_predictions_bonSigOnly.txt.gz";
			int rows=1;
			int cols=5;
			
			BufferedReader reader= FileUtils.createReader(fn);
			
			String line = null;
			int l = 0;
			while ((line=reader.readLine())!=null)
			{
				String[] eles = line.split("\t");
				String printBit = eles[0];
				if(cols<1)
					cols=eles.length;
				
				for(int c=1; c <  cols;c++)
				{
					printBit+="\t" + eles[c];
				}
				System.out.println(printBit);
				if(l>rows)
					break;
				l++;
			}
			System.out.println("Number of lines=\t" + l);
		}catch(Exception e){e.printStackTrace();}
	}

}
