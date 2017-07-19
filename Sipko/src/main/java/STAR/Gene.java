package STAR;

public class Gene 
{

	private String ensemblID;
	private String geneSymbol;
	
	public String chromosome;
	public int start;
	public int end;
	
	public String getEnsemblID()
	{
		return ensemblID;
	}

	public void setEnsemblID(String ensemblID)
	{
		this.ensemblID = ensemblID;
	}

	public String getGeneSymbol(boolean upperCase)
	{
		if(upperCase == true && this.geneSymbol!=null)
			return geneSymbol.toUpperCase();	
		return geneSymbol;
	}

	public void setGeneSymbol(String geneSymbol)
	{
		this.geneSymbol = geneSymbol;
	}

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
