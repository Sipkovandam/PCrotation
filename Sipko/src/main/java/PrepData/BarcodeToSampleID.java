package PrepData;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Hashtable;

import MatrixScripts.MyMatrix;
import Tools.FileUtils;

public class BarcodeToSampleID {
	//converts the fastq filename names in the counts file to the RNAseq id of interest
	//conversionFN uses Gerbens files 
	//static String conversionFN = "E:/Groningen/Data/RNAseq_clinic/Cardiomyopathy/1512_Cardio_RNAseq.txt";
	//static String countFN = "E:/Groningen/Test/JSON/ServerTest/Kallistocounts.txt.gz";
	
	static String conversionFN = "E:/Groningen/Data/RNAseq_clinic/RNAseq_Info.txt";
	static String countFN = "E:/Groningen/Data/RNAseq_clinic/5GPM/GRCh38/counts_GENES.txt";
		
	public static void main(String[] args) throws FileNotFoundException, IOException {
		checkArgs(args);
		Hashtable<String, String> conversion = getConversionHash();
		
		MyMatrix counts = new MyMatrix(countFN);
		for(int c = 0; c < counts.cols(); c++)
		{//160120_SN163_0694_AC8N1LACXX_L5_ATTCCT
			String[] eles = counts.colNames[c].split("_");
			int len = eles.length;
			if(len<5)
				continue;
			String identifier = eles[len-4]+"_"+eles[len-2]+"_"+eles[len-1];
			System.out.println("id = " + identifier);
			if(!conversion.containsKey(identifier))
				continue;
			counts.colNames[c] = conversion.get(identifier);
			System.out.println(identifier+ "\t" + counts.colNames[c]);
		}
		String outputFN = FileUtils.replaceEnd(countFN, "_RNAids.txt.gz");
		counts.write(outputFN);
		System.out.println(outputFN);
	}
	private static Hashtable<String, String> getConversionHash() throws FileNotFoundException, IOException {
		Hashtable<String, String> conversion = new Hashtable<String, String>();
		
		BufferedReader reader = FileUtils.createReader(conversionFN);
		String line = reader.readLine();//header line
		Hashtable<String, Integer> colIndex = FileUtils.makeHash(line);
		
		int barcode = colIndex.get("barcode");
		int lane = colIndex.get("lane");
		int run = colIndex.get("run");
		int name = colIndex.get("externalSampleID");
		while((line = reader.readLine())!=null)
		{
			String[] cells = line.split("\t");
			if(cells.length < barcode)
				continue;
			String runID = "";
			if(run< cells.length)
			{
				runID= cells[run];
				while(runID.length() < 4)
					runID="0"+runID;
			}
			System.out.println("key = " + runID+"_L"+cells[lane]+"_"+cells[barcode] + "\t" + cells[name]);
			conversion.put(runID+"_L"+cells[lane]+"_"+cells[barcode], cells[name]);
			//in some cases is different so just add that option as well:
			conversion.put("_"+cells[barcode]+"_L"+cells[lane], cells[name]);
		}
		return conversion;
	}
	static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko") && args.length < 1)
			return;
		if(args.length < 1)
		{
			System.out.println("Script requires the following argumetns:\n"
					+ "1. conversionFN=<conversionfn.txt> - File containing the barcodes (or whatever IDs are in the FastQfilenames) in 1st column and sampleID in 2nd\n"
					+ "2. fastqFolderName=<fastqFolderName.txt> - Folder where the corresponding FastQ files are\n");
			System.exit(1);
		}
		
		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase())
			{
//				var = new JSONutil<Vars>().read(var.JSON_FN, var);
				case "conversionfn":
					conversionFN =value;
					break;
				default:
					System.out.println("Incorrect argument supplied:\n"+ args[a] +"\nexiting");
					System.exit(1);
			}
		}
	}
}
