package Gdio;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.HashMap;

import Tools.FileUtils;
import Tools.Script;

public class QueryConverter extends Script<QueryConverter>
{
	String inputQueryFn = null;
	String queryWriteFn = null;
	String tableName_To_TableNameWithHashFn = null;
	
	
	@Override
	public void run()
	{
		try
		{
			HashMap<String,String> tableName_To_TableNameWithHash = FileUtils.readStringStringHash(tableName_To_TableNameWithHashFn);
			
			BufferedReader queryReader = FileUtils.createReader(inputQueryFn);
			String queryLine = null;
			StringBuilder newQuery = new StringBuilder();
			while((queryLine=queryReader.readLine())!=null)
			{
				queryLine=queryLine.replace(",", " ,").replace(")", " )").replace("(", "( ");
				String[] words=queryLine.split(" ");
				
				for(String word: words)
				{
					if(word.contains("."))
					{

						String newSubWord ="";
						String[] subWords = word.split("\\.");
						for(int w = 0; w < subWords.length; w++)
						{
							String subword= subWords[w];
							if(tableName_To_TableNameWithHash.containsKey(subword))
							{
								String convertedWord=tableName_To_TableNameWithHash.get(subword);
								newSubWord+=convertedWord;
							}
							else
								newSubWord+=subword;
							if(w < subWords.length-1)
								newSubWord+="\".\"";
						}
						word=newSubWord;
					}
					
					if(!word.matches("[A-Z]*|,|=|\\(|\\)"))
					{
						newQuery.append("\"");
					}
					
					if(tableName_To_TableNameWithHash.containsKey(word))
					{
						String newWord=tableName_To_TableNameWithHash.get(word);
						log(newWord);
						newWord=newWord.replace(".", "\".\"");
						newQuery.append(newWord);
					}
					else
						newQuery.append(word);
					
					if(!word.matches("[A-Z]*|,|=|\\(|\\)"))
					{
						newQuery.append("\"");
					}
					newQuery.append(" ");
				}
			}
			String result = newQuery.toString().replace(";\"", "\";");
			BufferedWriter resultWriter = FileUtils.createWriter(queryWriteFn);
			resultWriter.write(result+"\n");
			resultWriter.close();
			queryReader.close();
			log("New query is:\t" + result);
			log("Done query written to:\t" + queryWriteFn);
		}catch(Exception e){e.printStackTrace();}
		
		
	}

}

