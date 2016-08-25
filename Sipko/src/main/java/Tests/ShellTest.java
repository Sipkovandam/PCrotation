package Tests;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class ShellTest 
{

	public static void main (String[] args) throws IOException
	{
		String fileName = "E:/Groningen/Test/test.sh";
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(fileName)));
				
		writer.write("n=0\n"
				+ "while read line; do \n"
				+ "if [[ $line == \"Cronbach\"* ]]; then continue; fi\n"
				+ "compare=$(echo $line'<'0.7 | bc)\n"
				+ "if [[ compare -eq 1 ]]; then break; fi\n"
				+ "((n=$n+1))\n"
				+ "done < cronbachtest.txt\n"
				+ "echo $n\n");
		
		writer.close();
		System.out.println("done");
	}
}
