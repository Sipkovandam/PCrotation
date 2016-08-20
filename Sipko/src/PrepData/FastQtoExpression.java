package PrepData;

public class FastQtoExpression 
{
	//This class does the following:
	//1. Finds all fastQ files in a folder and subfolders
	//2. Determines which fastQ files are pairs
	//3. Maps the files using Kallisto
	//4. Creates a file describing the number of reads per sample mapped
	//5. Sums reads of transcritps to genes
	//6. Creates 1 large matrix for genes and one for transcripts
	//7. Creates 1 large matrix for samples >70% mapping
	
	
	public static class Var
	{
		String kallistoVersion = "";
		double mappingCutoff = 0.8;
	}
	
	public static void main(String args[])
	{
		
	}
	
}
