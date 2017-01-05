package TextEditing;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Hashtable;

import PCA.Matrix;
import Tools.FileUtils;

public class findAllGeneSymbols 
{

	public static void main (String[] args) throws FileNotFoundException, IOException
	{
		Matrix humanSymbols = new Matrix("D:/Sipko/06-10-2016/GeneNamesOnly.txt"); 
		
		BufferedReader reader = FileUtils.createReader("D:/Sipko/06-10-2016/ThesisText.txt");
		
		Hashtable<String,Integer> symbols = new Hashtable<String,Integer>();
		addSymbols(humanSymbols, symbols);
		
		String line = null;
		Hashtable<String,String> inThesis = new Hashtable<String,String>();
		while((line=reader.readLine())!=null)
		{
			String[] words = line.split(" |\t");
			for(int w=0; w<words.length;w++)
			{
				if(symbols.containsKey(words[w].toLowerCase()))
					inThesis.put(words[w],words[w]);
					
			}
		}
		Enumeration<String> wordsIn = inThesis.keys();
		
		int n = 0;
		while(wordsIn.hasMoreElements())
		{
			String word = wordsIn.nextElement();
			System.out.println(word);
			n++;
		}
		System.out.println(n);
	}
	
	

	private static void addSymbols(Matrix humanSymbols, Hashtable<String, Integer> symbols) {
		for(int r = 0; r < humanSymbols.rows(); r++)
		{
			symbols.put(humanSymbols.rowNames[r].toLowerCase(), 1);
		}
	}
}
