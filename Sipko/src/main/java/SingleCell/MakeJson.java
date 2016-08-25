package SingleCell;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class MakeJson {
	
	static String folderName = "";
	static String jsonName = "config.json";
	
	public static void main(String[] args) throws Exception
	{
		checkArgs(args);
		
		File files = new File(folderName);
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(jsonName)));
		for(File file : files.listFiles())
		{
			
//			writer.write("\"NUM_THREADS\": 32," + "\n");
//			writer.write("\"WINDOW\": [500, 5000]," + "\n");
//			writer.write("\"SOURCE_DIR\": "/path/to/source/"," + "\n");
//			writer.write("\"BASE_DIR\": "/path/to/fastq_files/"," + "\n");
//			writer.write("\"barcode_filenames\":" + "\n");
//			writer.write(barcodeFNs+ "\n");
//			writer.write(" \"read_filenames\":" + "\n");
//			writer.write("\"SAVE_DIR\": "/output/dir/for/save_data/"," + "\n");
//			writer.write("\"dmin\": 5," + "\n");
//			writer.write("\"BARCODE_LENGTH\": 14," + "\n");
//			writer.write("\"OUTPUT_DIR\": "/output/dir/for/singlecell_fastqs/"," + "\n");
//			writer.write("\"kallisto\":{" + "\n");
//			writer.write("	\"binary\": "/path/to/kallisto"," + "\n");
//			writer.write("	\"index\": "/path/to/kallisto_index.idx"," + "\n");
//			writer.write("	\"TCC_output\" : "/output/dir/for/TCC_output/"" + "\n");
//					
//					
//			writer.write("}");
//			writer.write("}");
		}
		writer.close();
	}
	static void checkArgs(String[] args) 
	{
		if(args.length < 1)
		{
			System.out.println("Wrong arguments");
			System.exit(1);
		}
		
		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase())
			{
				case "foldername":
					folderName =value;
					break;
				default:
					System.out.println("Incorrect argument supplied:\n"+ args[a] +"\nexiting");
					System.exit(1);
			}
		}
	}
}
