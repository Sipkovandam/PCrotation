package PCA;

import java.io.IOException;

public class MergeFiles 
{
	public static void main(String[] args) throws IOException
	{
		//Merges 2 files, keeping only the rows present in both files (alligns files based on rownames).
		
//		String inputFN = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/TESTexpression_Splits1_transposed.txt";
//		String additonFN = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/TESTexpression_Splits2_transposed.txt";
		
		String inputFN = "E:/Groningen/Data/bbmri/CountsGENES_0.7.txt";
		String additonFN = "E:/Groningen/Data/LifeLines/Phenotypes/CountsGENES_0.7.txt";
		
		boolean transposeInput = false;
		boolean transposeOutput = false;
				
		if(args.length!=0) 
		{
			inputFN = args[0];
			additonFN = args[1];
			if(args.length>2)
			{
				transposeInput = Boolean.parseBoolean(args[2]);
				transposeOutput = Boolean.parseBoolean(args[3]);
			}
		}
		
		MatrixStruct input = new MatrixStruct(inputFN);
		if(transposeInput)
			input.transpose();
		
		MatrixStruct addition = new MatrixStruct(additonFN);
		if(transposeOutput)
			addition.transpose();
		addition.keepRows(input);
		MatrixStruct result = input.mergeColumns(addition);
		
		if(transposeInput)
			result.transpose();
		
		result.write(inputFN.replace(".txt", "_merged.txt"));
	}

}
