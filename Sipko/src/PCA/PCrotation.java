package PCA;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.Arrays;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.TransformerException;

import no.uib.cipr.matrix.NotConvergedException;
import pca.PCA;

public class PCrotation 
{
	public static void main(String[] args) throws IOException, NotConvergedException, InterruptedException, ParserConfigurationException, TransformerException
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
//			case "rotatesamples":
//				RotateSample.main(argsToPass);
//				break;
//			case "correctforpcslude":
//				CorrectSampleForPCsLudeApproach.main(argsToPass);
//				break;
			case "getrows":
				GetRows.main(argsToPass);
				break;
			case "getcolumns":
				GetCols.main(argsToPass);
				break;
			case "getcols":
				GetCols.main(argsToPass);
				break;
			case "correctforpcs":
				PCcorrection.main(argsToPass);
				break;
			case "mergefiles":
				MergeFiles.main(argsToPass);
				break;
			case "transpose":
				Transpose.main(argsToPass);
				break;
			case "getoutliers":
				IdentifyDiseaseGenes.main(argsToPass);
				break;
			case "ranksamples":
				IdentifySimilarSamples.main(argsToPass);
				break;
			case "pcscores":
				PCA.scores(Paths.get(argsToPass[0]), Paths.get(argsToPass[1]));
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
				+ "1. geneEigenVectorsDirect"
				+ "2. correctForPCs\n"
				+ "3. getRows\n"
				+ "4. getColumns\n"
				+ "5. mergeFiles\n"
				+ "6. transpose\n"
				+ "7. getOutliers\n"
				+ "8. rankSamples\n"
				+ "Supply one of these arguments for futher information\n"
				+ "If you plan to use a large matrixes call this script like e.g.:"
				+ "java -jar -Xmx50g PCA.jar geneEigenVectors expressionFile.txt");
	}
	
}
