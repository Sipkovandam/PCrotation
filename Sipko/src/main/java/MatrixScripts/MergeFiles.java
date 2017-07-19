package MatrixScripts;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Set;

import Tools.FileUtils;
import Tools.Script;

public class MergeFiles extends Script<MergeFiles>
{
	//works a bit better with larger files, can be written even more efficiently if needed (memory wise)
	
	String fn1=null;
	String fn2=null;
	String writeFn=null;
	boolean keepAllFromFn1 = true;
	boolean keepAllFromFn2 = false;

	//String fn1="E:/Groningen/Data/PublicSamples/05-2016/Tessa/est_counts_nocancernocellline_Rlog_covariance_Chr1-22_add0.0BeforeGeoMean/CountsGENES_5GPM_46Samples/centered_transposed.txt";
	//String fn2="E:/Groningen/Data/PublicSamples/05-2016/Tessa/est_counts_nocancernocellline_Rlog_covariance_Chr1-22_add0.0BeforeGeoMean/CountsGENES_5GPM_46Samples/PC_1-300_DevidedBySTdevs.txt";
	
	@Override
	public void run() 
	{
		try
		{
			if(this.writeFn ==null)
				this.writeFn=FileUtils.addBeforeExtention(fn1, "_merged");
	
			HashMap<String, String> addRownamesToAddBit = FileUtils.readStringStringHash(fn2,true,1);
	
			BufferedReader f1Reader = FileUtils.createReader(fn1);
			BufferedReader f2Reader = FileUtils.createReader(fn2);
			String addHeader = f2Reader.readLine().split("\t",2)[1];
			BufferedWriter writer = FileUtils.createWriter(writeFn);
			String oldHeader = f1Reader.readLine();
			String newHeader = oldHeader+"\t"+addHeader;
			writer.write(newHeader+"\n");
			int newCols = newHeader.split("\t").length;
			int oldCols = oldHeader.split("\t").length;
			String placeHolder = "0";
			for(int n = oldCols; n < newCols-1 ; n++)
				placeHolder=placeHolder.concat("\t0");
			
			String line = null;

			Set<String> addedRows = addRownamesToAddBit.keySet();
			
			while((line = f1Reader.readLine())!=null)
			{
				String rowName = line.split("\t",2)[0];
				String addLine = addRownamesToAddBit.get(rowName);
				
				//only add once?
//				if(addedRows.contains(rowName))
//					addedRows.remove(rowName);
				
				if(addLine!= null)
					writer.write(line.concat("\t").concat(addLine).concat("\n"));
				else if (keepAllFromFn1)
					writer.write(line.concat("\t").concat(placeHolder).concat("\n"));
			}
			
			
			
			if(keepAllFromFn2)
			{
				String placeHolder2="";
				for(int n = 0; n < oldCols-1; n++)
					placeHolder2=placeHolder2.concat("\t0");

				for(String geneName: addedRows)
				{
					String addLine = addRownamesToAddBit.get(geneName);
					writer.write(geneName.concat(placeHolder2).concat("\t").concat(addLine).concat("\n"));
				}
			}
			writer.close();
			f1Reader.close();
			f2Reader.close();
		}catch(Exception e){e.printStackTrace();}
	}

	public String getFn1()
	{
		return fn1;
	}

	public void setFn1(String fn1)
	{
		this.fn1 = fn1;
	}

	public String getFn2()
	{
		return fn2;
	}

	public void setFn2(String fn2)
	{
		this.fn2 = fn2;
	}

	public String getWriteFn()
	{
		return writeFn;
	}

	public void setWriteFn(String writeFn)
	{
		this.writeFn = writeFn;
	}

	public boolean isKeepAllFromFn1()
	{
		return keepAllFromFn1;
	}

	public void setKeepAllFromFn1(boolean keepAllFromFn1)
	{
		this.keepAllFromFn1 = keepAllFromFn1;
	}

	public boolean isKeepAllFromFn2()
	{
		return keepAllFromFn2;
	}

	public void setKeepAllFromFn2(boolean keepAllFromFn2)
	{
		this.keepAllFromFn2 = keepAllFromFn2;
	}

}
