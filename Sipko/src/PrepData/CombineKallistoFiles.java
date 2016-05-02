package PrepData;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Hashtable;

import PCA.Matrix;
import pca.PCA;

public class CombineKallistoFiles 
{
	static int column = 2; //2 is read counts, 3 is TPM values in kallisto output
	
	public static void main (String[] args)
	{
//		String folder = "E:/Groningen/Data/LifeLines/Phenotypes/Kallisto mapping counts/";//Takes 10 minutes to combine these 1200 files
//		String folderMappingPercentage = "E:/Groningen/Data/LifeLines/quantification/";
		
		String folderFN = "E:/Groningen/Data/bbmri/kallisto/";//Takes 10 minutes to combine these 1200 files
		String folderMappingPercentage = "E:/Groningen/Data/bbmri/quantification/";
		//String folder = "E:/Groningen/Data/5GPM/kallisto/";//Takes 10 minutes to combine these 1200 files
		//String folderMappingPercentage = "E:/Groningen/Data/5GPM/kallisto/jobs_QT/";
		double threshold = 0.7;
		boolean sensible = true;//if using niek/freerks file output structure this should be false (see what I did there? :D)
		String transcriptsToGenesFile = "E:/Groningen/Data/Annotation/hg19.v75.cdna.all.enst2ensg.txt";
		
		//String dir = System.getProperty("user.dir");
		if(!System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
		{	
			if(args.length < 4)
			{
				printMsg();
			}
			else
			{
				folderFN = args[0];
				folderMappingPercentage = args[1];
				threshold = Double.parseDouble(args[2]);
				transcriptsToGenesFile = args[3];
				column = Integer.parseInt(args[4]);
			}
		}
		
		File file = new File(folderFN);
		String parentPath = file.getAbsoluteFile().getParent();
		String writeFileName = parentPath+"/Counts.txt";
		if(column == 3)
			writeFileName = writeFileName.replace(".txt", "_TPM_.txt");
		//Combining files
		System.out.println("Combining expression files in folder: " + folderFN);
		File fldr = new File(folderFN);
		if(!fldr.exists())
		{
			System.out.println("This folder does not exist: " + folderFN);
			System.exit(1);
		}
		File[] fileNames = fldr.listFiles();		
		Matrix output = null;
		//check for missing files first...
		int missing = 0;
		BufferedWriter fw;
		int nFiles = fileNames.length;
		try 
		{
			fw = new BufferedWriter(new FileWriter(new File(parentPath+"/missing.txt")));
			for(int f = 0; f < nFiles; f++)
			{
				String fileName = fileNames[f].getAbsolutePath()+"/abundance.tsv";
				File fl = new File(fileName);
				if(!fl.exists())
				{
					System.out.println("The following file is missing and excluded: " + fileName);
					missing++;
					fw.write(fileName+"\n");
				}
			}
			fw.close();
		} catch (IOException e) {e.printStackTrace();}
		
		if(missing!=0) 
			System.out.println("A total of:" + missing + " files is missing");
		
		int outC=0;
		for(int f = 0; f < nFiles; f++)
		{
			String fileName = fileNames[f].getAbsolutePath()+"/abundance.tsv";
			File fl = new File(fileName);
			if(!fl.exists())
				continue;
			
			if(f%100==0)
			{
				System.out.println("File: " + f + "/" + fileNames.length);
				System.out.println(fileName);
				PCA.log(Double.toString(((double) f)/((double)fileNames.length)));
			}
			
			Matrix counts = new Matrix(fileName);
			
			if(output == null)
			{
				output = new Matrix(counts.rowNames.length, nFiles-missing);
				output.rowNames = counts.rowNames;
			}
			File tempName = new File(fileName);
			String[] split = tempName.getParent().split("\\\\");
			output.colNames[outC] =  split[split.length-1];
			for(int x = 0; x < counts.rowNames.length;x++)
			{
				if(output.rowNames[x].compareTo(counts.rowNames[x]) != 0)
				{
					System.out.println("RowNames are diffferent, need to use a index hash!");//just use Matrix.removeRows(Matrix a) function (removes rows not present in both files and alligns the others)
					System.out.println("output.rowNames[x] = " + output.rowNames[x] + " counts.rowNames[x] = " + counts.rowNames[x] + " x = " + x);
				}
				//counts.print(4,-1);
				//System.out.println("X =" + x + " " +counts.values[x].length);
				output.values[x][outC] = counts.values[x][column];
			}
			outC++;
		}
		output.print(3, 5);
		System.out.println("WriteFN= " + writeFileName);
		output.write(writeFileName, -1);
		System.out.println("Files combined into 1 expression file: " + writeFileName);

		//sumTotranscripts
		System.out.println("Summing transcripts to genes");
		String[] arguments = new String[]{writeFileName,transcriptsToGenesFile};
		Matrix geneExpressionCounts = SumTranscriptsToGenes.run(writeFileName,transcriptsToGenesFile);//(also saves a file)
		System.out.println("Done summing transcripts to genes");

		//Create a file that reports all the percentages of mapping reads
		System.out.println("Creating mapping % file");
		Matrix percentages = null;
		
		if(!folderMappingPercentage.toLowerCase().contains("noquality"))
		{
			if(!sensible)
				percentages = getPercentagesNiekFreerk(folderMappingPercentage, folderFN);
			if(sensible)
				percentages = getPercentages(folderMappingPercentage, "mappingPercentages.txt");
			
			String wn = parentPath+"/mappingPercentages.txt";
			percentages.print(5,5);
			percentages.write(wn);
			System.out.println("Percentages mapping file created at: " + wn);
		}
		else//if no quality files are available just combine all files
		{
			File folder = new File(folderFN);
			percentages = new Matrix(folder.listFiles().length,1);
			for(int r = 0; r < percentages.rowNames.length; r++)
				percentages.values[r][0] = 1;
		}
		
		
		
		

		DecimalFormat formatter = new DecimalFormat("#.##");
		String thresholdSaveName = parentPath+"/CountsGENES_"+formatter.format(threshold)+".txt";
		if(column == 3)
			thresholdSaveName = thresholdSaveName.replace(".txt", "_TPM_.txt");
		//remove all samples under a certain mapping threshold.
		System.out.println("Removing samples under threshold: " + threshold);
		int n = 0;
		for(int f = 0; f < percentages.rowNames.length ;f++)
		{
			if(percentages.values[f][0] > threshold && percentages.values[f][0] >0)
				n++;
		}
		
		Matrix thresholdExpressionSamples = new Matrix(geneExpressionCounts.rowNames.length,n);
		//change percentages matrix names to have the same names as the colNames of the files
		for(int r = 0; r < percentages.rowNames.length; r++)
		{
			if(!percentages.rowNames[r].contains("abundance.tsv"))
				continue;
			File tempName = new File(percentages.rowNames[r].replace("abundance.tsv", ""));
			percentages.rowNames[r] = tempName.getName();
		}
		
		Hashtable<String, Integer> colNamesToNumbers = new Hashtable<String,Integer>();
		String folderName2 = new File(geneExpressionCounts.colNames[1]).getParent();
		for(int c = 0; c < geneExpressionCounts.colNames.length; c++)
		{
			File tempName = new File(geneExpressionCounts.colNames[c]);
			colNamesToNumbers.put(tempName.getName().substring(0,tempName.getName().replace("_internal_id", "").lastIndexOf('_')),c);
		}
		
		int missing2 = 0;
		for(int r = 0; r < geneExpressionCounts.rowNames.length ;r++)
		{
			int col = 0;
			for(int c = 0; c < percentages.rowNames.length; c++)
			{
				if(percentages.values[c][0] > threshold)
				{
					if(colNamesToNumbers.get(percentages.rowNames[c]) == null)
					{
						if(r == 0)
						{
							System.out.println("Warning, the counts for the following sample was not found (and is thus exlcuded): " + percentages.rowNames[c]);// + " hash = " + colNamesToNumbers.keys().nextElement());
							missing2++;
						}
						continue;
					}
					thresholdExpressionSamples.colNames[col] = percentages.rowNames[c];
					thresholdExpressionSamples.values[r][col] = geneExpressionCounts.values[r][colNamesToNumbers.get(percentages.rowNames[c])];
					col++;
				}
			}
		}
		thresholdExpressionSamples.rowNames = geneExpressionCounts.rowNames;
		System.out.println("Number of missing samples= " + missing2 + " Number of samples for which mapping counts are unknown: " + missing);
		
		thresholdExpressionSamples.write(thresholdSaveName,-1);
		System.out.println("Samples removed under: " + threshold + " file written to : " +thresholdSaveName);
		
	}

	private static Matrix getPercentages(String folderMappingPercentage, String percentFN) {
		
		try {
			SearchFilesInDirectories.main(new String[]{folderMappingPercentage+"/",percentFN});
		} catch (IOException e) {e.printStackTrace();}
		String percentagesFN = folderMappingPercentage+"/"+ percentFN+".txt";
		Matrix percentagesFNs = new Matrix(percentagesFN);
		Matrix output = new Matrix(percentagesFNs.rowNames.length,1);
		output.rowNames = percentagesFNs.rowNames;
		output.colNames[0] = "Mapping percentages";
		for(int r = 0; r < percentagesFNs.rowNames.length; r++)
		{
			String fileName = percentagesFNs.rowNames[r];
			String[] toGet = new String[]{"[quant] processed "//gets line with the number of mapping reads
										 ,"[quant] will process"};//gets the line with the fileName
			String lines[] = getLines(toGet,fileName);
			if(lines[0] == null)
			{
				System.out.println("This file does not have any information on the number of reads mapped:" + fileName);
				continue;
			}
			String[] eles = lines[0].split(" ");
			double mapping = Double.parseDouble(eles[4].replace(",", ""));
			double total = Double.parseDouble(eles[2].replace(",", ""));
			double percent = mapping/total;
			output.values[r][0] = percent;
		}
		
		return output;
	}

	private static Matrix getPercentagesNiekFreerk(String folderMappingPercentage, String folderReadCounts) 
	{
		//Create a hashtable that has the readcount FileName based on the input file for the mapping (the .fq.gz file)
		Hashtable<String, String> inputToCountFilename = new Hashtable<String,String>();
		File folder = new File(folderReadCounts);
		File[] files = folder.listFiles();
		
		//make hash that converts 
		for(int f = 0; f < files.length; f++)
		{
			if(!files[f].isDirectory())
				continue;
			String folderName = files[f].getAbsolutePath();
			String fileName = files[f].getName().substring(0,files[f].getName().replace("_internal_id", "").lastIndexOf('_'));
			//System.out.println("FN =" + fileName);
			inputToCountFilename.put(fileName, folderName+"/abundance.tsv");
		}
		
		File folderPercentage = new File(folderMappingPercentage);
		
		files = folderPercentage.listFiles();
		int len = 0;
		for(int f = 0; f < files.length; f++)
		{
			String fileName = files[f].getAbsolutePath();
			if(!fileName.endsWith(".err"))
				continue;
			len++;
		}
		System.out.println("len = " + len);
		Matrix output = new Matrix(len,1);
		
		output.colNames[0] = "Percentage";
		int out = 0;
		for(int f = 0; f < files.length; f++)
		{
			String fileName = files[f].getAbsolutePath();
			if(!fileName.endsWith(".err"))
				continue;
			String[] toGet = new String[]{"[quant] processed "//gets line with the number of mapping reads
										 ,"[quant] will process"};//gets the line with the fileName
			String lines[] = getLines(toGet,fileName);
			if(lines[0] == null)
			{
				System.out.println("This file does not have any information on the number of reads mapped:" + fileName);
				continue;
			}
			String[] eles = lines[0].split(" ");
			double mapping = Double.parseDouble(eles[4].replace(",", ""));
			double total = Double.parseDouble(eles[2].replace(",", ""));
			double percent = mapping/total;
			File tempFileName = new File(lines[1]);
			String lineName = tempFileName.getName();
			String name = lineName.substring(0,lineName.replace("_internal_id", "").lastIndexOf('_'));
			String targetName = inputToCountFilename.get(name);
			if(targetName == null)
			{
				System.out.println("name =" + name);
				System.out.println("size nameHash = " + inputToCountFilename.size());
				System.out.println("nameHash =" + inputToCountFilename.keys().nextElement());
			}
			output.rowNames[out] = name;
			output.values[out][0] = percent;
			out++;
		}
		return output;
	}

	private static String[] getLines(String[] toGet, String fileName) 
	{
		String[] output = new String[toGet.length];
		try 
		{
			BufferedReader reader = new BufferedReader(new FileReader(new File(fileName)));
			String line = null;
			while((line = reader.readLine())!=null)
			{
				for(int e = 0; e < toGet.length; e++)
				{
					if(line.contains(toGet[e]))
					{
						output[e] = line;
					}
				}
			}
			
			reader.close();
		} catch (FileNotFoundException e) {e.printStackTrace();} catch (IOException e) {e.printStackTrace();}

		return output;
	}

	private static void printMsg() {
		System.out.println("Function takes 1 argumentes\n"
				+ "1. Foldername where the kallisto count files are located"
				+ "2. noquality, if no quality files are available");
		
	}
}
