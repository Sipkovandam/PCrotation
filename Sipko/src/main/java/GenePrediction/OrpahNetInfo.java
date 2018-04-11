package GenePrediction;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

import Tools.FileUtils;

public class OrpahNetInfo
{
	HashMap<String,HashSet<String>> disease_To_Genes =  new HashMap<String,HashSet<String>>();
	HashMap<String,HashSet<String>> disease_To_HpoTerms =  new HashMap<String,HashSet<String>>();
	HashMap<String, HashSet<String>> genes_To_diseases = new HashMap<String, HashSet<String>>();
	HashMap<String, HashSet<String>> genes_To_HpoTerms = new HashMap<String, HashSet<String>>();
	
	OrpahNetInfo(String fn) throws FileNotFoundException, IOException
	{
		readOphaNetFile(fn);
	}
	
	public void readOphaNetFile(String fn) throws FileNotFoundException, IOException
	{
		 BufferedReader orphaReader = FileUtils.createReader(fn);
		 
		 String line = orphaReader.readLine();
		 
		 while((line=orphaReader.readLine())!=null)
		 {
			 String eles[] = line.split("\t");
			 String orphaDisease = eles[0];
			 String geneSymbol = eles[1];
			 String geneId  = eles[2];
			 String hpoTermId = eles[3];
			 String hpoTermName = eles[4];
			 
			 
			 //add element (first argument) to hash (second argument)
			 addElementToHashMap(geneId,orphaDisease, disease_To_Genes);
			 addElementToHashMap(hpoTermId, orphaDisease, disease_To_HpoTerms);
			 addElementToHashMap(orphaDisease, geneId, genes_To_diseases);
			 addElementToHashMap(hpoTermId, geneId, genes_To_HpoTerms);
			 
		 }
	}

	private void addElementToHashMap(	String hashName,
									String elementToAdd, HashMap<String,HashSet<String>> hashMap)
	{
		HashSet<String> elements = hashMap.get(hashName);
		if(elements==null)
			elements = new HashSet<String>();
		elements.add(elementToAdd);
		hashMap.put(hashName, elements);
	}

}
