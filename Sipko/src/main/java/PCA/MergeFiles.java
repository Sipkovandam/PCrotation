package PCA;

import java.io.IOException;

public class MergeFiles
{
	public static void main(String[] args) throws IOException
	{
		//Merges 2 files, keeping only the rows present in both files (alligns files based on rownames).

		//String inputFN = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/TESTexpression_Splits1_transposed.txt";
		//String additonFN = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/TESTexpression_Splits2_transposed.txt";

		String inputFN = "E:/Sipko/Thesis resubmission/FDR analysis/TopXpvalues_genes.out";
		String additonFN = "E:/Sipko/Thesis resubmission/Reactome Files/EnsemblGeneToType.txt";

		boolean transposeInput = false;
		boolean transposeOutput = false;
		boolean keepAll = true;

		if (args.length != 0)
		{
			inputFN = args[0];
			additonFN = args[1];
			if (args.length > 2)
			{
				transposeInput = Boolean.parseBoolean(args[2]);
				transposeOutput = Boolean.parseBoolean(args[3]);
			}
		}

		MatrixStruct input = new MatrixStruct(inputFN);
		if (transposeInput)
			input.transpose();

		MatrixStruct addition = new MatrixStruct(additonFN);
		if (transposeOutput)
			addition.transpose();
		if (keepAll)
			addition.keepRows1Matrix(input);
		else
			addition.keepRows(input);
		MatrixStruct result = input.mergeColumns(	addition,
													keepAll);

		if (transposeInput)
			result.transpose();

		result.write(inputFN.replace(	".txt",
										"_merged.txt"));
	}

}
