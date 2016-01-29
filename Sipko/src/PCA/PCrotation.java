package PCA;

import java.io.IOException;
import java.util.Arrays;

import no.uib.cipr.matrix.NotConvergedException;

public class PCrotation 
{
	public static void main(String[] args) throws IOException, NotConvergedException
	{
		if(args.length == 0)
		{
			printUsage();
		    System.exit(1);
		}
		String[] argsToPass = Arrays.copyOfRange(args,1,args.length);
		switch (args[0].toLowerCase()) 
		{
			case "geneeigenvectors":
				CreateGeneEigenvectorFile.main(argsToPass);
				break;
			case "rotatesamples":
				RotateSample.main(argsToPass);
				break;
			case "correctforpcs":
				CorrectForPCs.main(argsToPass);
				break;
			case "correctforpcslude":
				CorrectSampleForPCsLudeApproach.main(argsToPass);
				break;
			case "getrows":
				GetRows.main(argsToPass);
				break;
			default:
				printUsage();
			    System.exit(1);
		}
	}

	private static void printUsage() 
	{
		System.out.println("This script can be called with the following arguments:\n"
				+ "1. geneEigenVectors\n"
				+ "2. rotateSamples\n"
				+ "3. correctForPCs\n"
				+ "Supply one of these arguments for futher information\n"
				+ "If you plan to use a large matrixes call this script like e.g.:"
				+ "java -jar -Xmx100g PCA.jar geneEigenVectors expressionFile.txt");
	}
	
}
