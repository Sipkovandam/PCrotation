package PrepData;

import java.io.IOException;

import PCA.MatrixStruct;

public class CorrelationToPopulateDBlayout 
{
	//converts the correlation matrix layout so it can be used for: populateCoregulationDBTXT.js
	//this mainly entails sorting the correlation matrix gene names in the same way as the <format.txt> file
	//additionally it removes all genes not present in <format.txt>
	
	public static void main(String[] args) throws IOException
	{
	String correlationFN = "E:/Groningen/Data/Juha/node_modules/gene_correlation.txt";
	String keepGenesFN = "E:/Groningen/Data/Juha/node_modules/GeneNameHGNCV71_sortedGeneNames.txt";
	String writeFN = null;
	
	if(args.length <2)
		checkArgs(args);
	for(int a = 0; a < args.length; a++)
	{
		String arg = args[a].split("=")[0];
		String value = args[a].split("=")[1];
		switch (arg.toLowerCase()){
			case "filename":
				correlationFN =value;
				break;
			case "writefn":
				writeFN = value;
				break;
			case "keepgenesfn":
				keepGenesFN = value;
				break;
			default:
				checkArgs(args);
				System.out.println("Incorrect argument supplied; exiting");
				System.exit(1);
		}
	}
	if(writeFN == null)
		writeFN = correlationFN.replace(".txt", "_correctOrder.txt");
	
	MatrixStruct correlation = new MatrixStruct(correlationFN);
	MatrixStruct keepGenes = new MatrixStruct(keepGenesFN);
	adjust(correlation, keepGenes, writeFN);
}

public static void adjust(MatrixStruct correlation, MatrixStruct keepGenes, String writeFN) throws IOException 
{
	keepGenes.keepRows(correlation);
	correlation.transpose();
	keepGenes.keepRows(correlation);
	correlation.write(writeFN);
	System.out.println("File written to:" + writeFN);
}
static void checkArgs(String[] args) 
{
	if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
		return;
	System.out.println("Wrong arguments, this script takes the following arguments:\n"
			+ "1. filename=<path/correlation.txt> - filename of the correlation matrix over the genes for the PCs > 0.7 cronbach alpha\n"
			+ "2. keepGenesFN=<path/GeneNamesOrder.txt> - Data file (must not contain Strings in any other fields then header or rowname fields)\n"
			+ " containing the genes to keep in the row names (assumes a header row as well) \n"
			+ "3. writeFN=<path/writeFN.txt> - filename of the file to write");
	System.exit(1);
}
}
