package Tests;

import Tools.Run;

public class RunTest
{
	
	public static void main(String[] args) throws Exception
	{
		//args = new String[] { "E:/Groningen/Test/JSON/PCApipeline/config_PCA.PCApipeline_2017-02-10.json" };
		
		String testClass = "GeneNetwork.WebsiteMatrixCreator";
		args = new String[]{"getjson",testClass,"E:/Groningen/Test/"+testClass+"/"};
		 
//		args = new String[]{"E:/Groningen/Test/GeneNetwork.ChildHpoBasedZscorePredictor/config_GeneNetwork.ChildHpoBasedZscorePredictor_2018-03-12.json"};		
//		args = new String[]{"E:/Groningen/Test/GeneNetwork.SubHpoPredictor/config_GeneNetwork.SubHpoPredictor_2018-03-09.json"};
		
//		args = new String[]{"E:/Groningen/Test/GenePrediction.GenePredictor/config_GenePrediction.GenePredictor_2018-04-04.json"};
		
		//String testClass = "GeneNetwork.CovariateAdjuster";
		//args = new String[]{"getjson",testClass,"E:/Groningen/Test/"+testClass+"/"};
			
//		String testClass = "GeneNetwork.GetGeneRanks";
//		args = new String[]{"E:/Groningen/Test/GeneNetwork.GetGeneRanks/BenchmarkSamples3/config_GeneNetwork.GetGeneRanks_2018-03-06.json"};
		
//		args = new String[]{"E:/Groningen/Test/MatrixScripts.Transpose/config_MatrixScripts.Transpose_2017-06-05.json"};
//		args = new String[]{"E:/Groningen/Test/Kallisto._Kallisto_Pipeline/config_Kallisto._Kallisto_Pipeline_2017-10-19.json"};
		//args = new String[]{"E:/Groningen/Test/GeneNetwork.WebsiteMatrixCreator/config_GeneNetwork.WebsiteMatrixCreator_GO_C.json"};
		args = new String[]{"E:/Groningen/Data/Annotation/GRCh38/PathwaysPatrick/config_MatrixScripts.GetRows_GO_F.json"};
		
		Run.run(args);
	}
}
