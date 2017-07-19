package PizzlyClasses;

import java.util.HashMap;

import Tools.Script;

public class Transcript extends Script<Transcript>
{
	String fasta_record =null;
	HashMap<String,String> transcriptA = null;
	HashMap<String,String> transcriptB = null;
	int support = 0;
	int[] reads = null;
}
