package Tools;

import java.io.BufferedReader;
import java.util.HashMap;

public class SampleSheetGafParser extends Script<SampleSheetGafParser>
{
	/**
	 * 
	 */
	private static final long serialVersionUID = -6648836498592663166L;
	public transient HashMap<String,HashMap<String,String>> sampleSheet = new HashMap<String,HashMap<String,String>>();//<columnName,HashMap<sampleName,entry>
	private String sampleSheetFn = null;
	
	public void run()
	{
		try
		{
			BufferedReader sampleSheetReader = FileUtils.createReader(sampleSheetFn);
			HashMap<Integer,String> colNumberToColname= new HashMap<Integer, String>();
			HashMap<String,Integer> colNameToColNumber= new HashMap<String, Integer>();
			
			String header=sampleSheetReader.readLine();
			p("header = " +header);
			String[] colNames = header.split(separator);
			for(int c = 0; c < colNames.length; c++)
			{
				String colName = colNames[c];
				sampleSheet.put(colName, new HashMap<String,String>());
				colNumberToColname.put(c, colName);
				colNameToColNumber.put(colName, c);
			}
			
			String line = null;
			while((line=sampleSheetReader.readLine())!=null)
			{
				String[] dataPerCol=line.split(separator);
				int sequencingStartDateCol = colNameToColNumber.get("sequencingStartDate");
				int sequencerCol = colNameToColNumber.get("sequencer");
				int runCol = colNameToColNumber.get("run");
				int flowcellCol = colNameToColNumber.get("flowcell");
				int lane = colNameToColNumber.get("lane");
				int barcode = colNameToColNumber.get("barcode");
				
				String sampleName = dataPerCol[sequencingStartDateCol]+"_"+
						dataPerCol[sequencerCol]+"_"+
						addZeros(dataPerCol[runCol])+"_"+
						dataPerCol[flowcellCol]+"_L"+
						dataPerCol[lane]+"_"+
						dataPerCol[barcode];
				
				for(int c = 0; c < colNames.length; c++)
				{
					String colName = colNumberToColname.get(c);
					HashMap<String,String> colEntries = sampleSheet.get(colName);
					colEntries.put(sampleName, dataPerCol[c]);
				}
			}
			
		}catch(Exception e)
		{
			e.printStackTrace();
		}
	}
	public void loadSampleSheet()
	{
		this.run();
	}

	private String addZeros(String run)
	{
		while(run.length()<3)
			run = "0"+run;
		return run;
	}
	public String getSampleSheetFn()
	{
		return sampleSheetFn;
	}
	public void setSampleSheetFn(String sampleSheetFn)
	{
		this.sampleSheetFn = sampleSheetFn;
	}
	public String getSeparator()
	{
		return separator;
	}
	public void setSeparator(String separator)
	{
		this.separator = separator;
	}

	private String separator=",";
}
