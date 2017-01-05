package STAR;

public class Gene 
{

	public String ensemblID;
	public String geneSymbol;
	
	public String chromosome;
	public int start;
	public int end;
	
	Gene(String ensemblID, String geneSymbol, String chromosome, int start, int end)
	{
		this.ensemblID = ensemblID;
		this.geneSymbol = geneSymbol;
		this.chromosome = chromosome;
		this.start = start;
		this.end = end;
	}

	public boolean inGene(int position) {
		if(position>start && position < end)
			return true;
		return false;
	}
	
}
