package Tools;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.Arrays;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.TransformerException;

import Analyses.AveragesPerRow;
import Analyses.CorrelateFiles;
import JuhaPCA.PCA;
import Kallisto.CombineKallisto;
import Kallisto.FastQtoExpression;
import Kallisto.Slurm;
import Kallisto.KeepThresholdSamples;
import Kallisto.RemoveBadSamples;
import Kallisto.SumTranscriptsToGenes;
import PCA.PCApipeline;
import PCA.GetCols;
import PCA.GetRows;
import PCA.IdentifyDiseaseGenes;
import PCA.IdentifySimilarSamples;
import PCA.LogTransform;
import PCA.MergeFiles;
import PCA.RLog;
import PCA.RlogLargeMatrix;
import PCA.RlogLargeMatrix_Main;
import PCA.Transpose;
import PCA.Zscore;
import PrepData.GetSamplesWithEmptyCells;
import STAR.SpliceMerger;
import STAR._STAR_Pipeline;
import no.uib.cipr.matrix.NotConvergedException;

public class Toolkit
{
	//Script that refers to all the other scripts. This can be called...

	public static void main(String[] args) throws Exception
	{
		if (args.length == 0)
		{
			printUsage();
			System.exit(1);
		}
		String[] argsToPass = Arrays.copyOfRange(	args,
													1,
													args.length);
		switch (args[0].toLowerCase())
		{
		case "json":
			runScript(argsToPass);
			break;
		case "run":
			Run.run(argsToPass);
			break;
		case "getsampleswithemptycells":
			GetSamplesWithEmptyCells.main(argsToPass);
			break;
		case "rlog":
			RLog.main(argsToPass);
			break;
		case "rloglarge":
			RlogLargeMatrix_Main.main(argsToPass);
			break;
		case "getrows":
			GetRows.main(argsToPass);
			break;
		case "getcolumns":
			GetCols.main(argsToPass);
			break;
		case "getcols":
			GetCols.main(argsToPass);
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
		//			case "pcscores":
		//				PCA.scores(Paths.get(argsToPass[0]), Paths.get(argsToPass[1]));
		//				break;
		case "removebadsamples":
			RemoveBadSamples.main(argsToPass);
			break;
		case "combinekallisto":
			CombineKallisto.main(argsToPass);
			break;
		case "kallistoslurm":
			new Slurm().run(argsToPass);
			break;
		case "keepthresholdsamples":
			KeepThresholdSamples.main(argsToPass);
			break;
		case "fastqtoexpression":
			FastQtoExpression.main(argsToPass);
			break;
		case "matrixstats":
			MatrixStats.main(argsToPass);
			break;
		case "sumtranscriptstogenes":
			SumTranscriptsToGenes.main(argsToPass);
			break;
		case "averagesperrow":
			AveragesPerRow.main(argsToPass);
			break;
		case "zscores":
			Zscore.main(argsToPass);
			break;
		case "correlatefiles":
			CorrelateFiles.main(argsToPass);
			break;
		case "log2":
			LogTransform.main(argsToPass);
			break;
		default:
			printUsage();
			System.exit(1);
		}

	}

	private static void runScript(String[] argsToPass) throws IOException
	{
		String jsonName = argsToPass[0];
		String scriptName = FileUtils.getLine(	jsonName,
												"\"className\": \"")
				.split("\"")[3];

		switch (scriptName)
		{
		case "RlogLargeMatrix":
			RlogLargeMatrix rlogLargeMatrix = new JSONutil<RlogLargeMatrix>().read(	jsonName,
																					new RlogLargeMatrix());
			rlogLargeMatrix.run();
			break;

		}
	}

	private static void printUsage()
	{
		System.out.println("This script can be called with the following arguments:\n" + "1.  PCAsteps\n" + "2.  PCcorrection\n" + "3.  getRows\n" + "4.  getColumns\n" + "5.  mergeFiles\n" + "6.  transpose\n" + "7.  getOutliers\n" + "8.  rankSamples\n" + "9.  rLog\n" + "10. rLogLarge\n" + "11. matrixStats\n" + "12. SumTranscriptsToGenes" + "13. GeneToTranscript\n" + "14. AveragesPerRow\n" + "15. zscores\n" + "16. correlatefiles\n" + "17. correlationvsexpressionplots\n" + "18. log2\n" + "19. fastqmappingstar\n" + "20. SpliceSites\n" + "Supply one of these arguments for futher information\n" + "If you plan to use a large matrixes call this script like e.g.:" + "java -jar -Xmx50g PCA.jar geneEigenVectors expressionFile.txt");
	}

}
