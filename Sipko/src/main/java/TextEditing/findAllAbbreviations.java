package TextEditing;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Enumeration;
import java.util.Hashtable;

import MatrixScripts.MyMatrix;
import Tools.FileUtils;

public class findAllAbbreviations {

	public static void main(String[] args) throws FileNotFoundException, IOException 
	{
		MyMatrix already = new MyMatrix("D:/Sipko/06-10-2016/AlreadyInAbreviations.txt");
		
		BufferedReader reader = FileUtils.createReader("D:/Sipko/06-10-2016/ThesisText.txt");
		
		String line = null;
		Hashtable<String,String> inThesis = new Hashtable<String,String>();
		while((line=reader.readLine())!=null)
		{
			String[] words = line.split(" |\t");
			for(int w=0; w < words.length;w++)
			{
				if(already.getRowHash().containsKey(words[w]) 
						|| words[w].contains("~") 
						|| words[w].contains("|")
						|| words[w].contains("$") 
						|| words[w].contains("0") 
						|| words[w].length() > 10)
					continue;
				int capitals = countCapitals(words[w]);
				if(capitals>1)
					inThesis.put(words[w].replace(",", "").replace(".", ""), words[w].replace(",", "").replace(".", ""));
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

	private static int countCapitals(String string) {
		int upperCase = 0;
		for (int k = 0; k < string.length(); k++) {
		    /**
		     * The methods isUpperCase(char ch) and isLowerCase(char ch) of the Character
		     * class are static so we use the Class.method() format; the charAt(int index)
		     * method of the String class is an instance method, so the instance, which,
		     * in this case, is the variable `input`, needs to be used to call the method.
		     **/
		    // Check for uppercase letters.
		    if (Character.isUpperCase(string.charAt(k))) upperCase++;

//		    // Check for lowercase letters.
//		    if (Character.isLowerCase(string.charAt(k))) lowerCase++;
		}
		return upperCase;
	}

}
