package PizzlyClasses;

import java.util.HashMap;
import Tools.Script;

public class Fusion extends Script<Fusion>
{
	HashMap<String,String> geneA = null;
	HashMap<String,String> geneB = null;
	
	int paircount = 0;
	int splitcount = 0;
	
	Transcript[] transcripts = null;
	
	Readpair[] readpairs = null;

	public HashMap<String, String> getGeneA()
	{
		return geneA;
	}

	public void setGeneA(HashMap<String, String> geneA)
	{
		this.geneA = geneA;
	}

	public HashMap<String, String> getGeneB()
	{
		return geneB;
	}

	public void setGeneB(HashMap<String, String> geneB)
	{
		this.geneB = geneB;
	}

	public int getPaircount()
	{
		return paircount;
	}

	public void setPaircount(int paircount)
	{
		this.paircount = paircount;
	}
}