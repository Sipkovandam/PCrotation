package Tests;

import Tools.Run;

public class RunTest
{
	
	public static void main(String[] args) throws Exception
	{
		//args = new String[] { "E:/Groningen/Test/JSON/PCApipeline/config_PCA.PCApipeline_2017-02-10.json" };
		
		String testClass = "Kallisto.SumTranscriptsToGenes";
		args = new String[]{"getjson",testClass,"E:/Groningen/Test/"+testClass+"/"};
		 
//		args = new String[]{"E:/Groningen/Test/Analyses.OutdatedHpoGenesRetriever/config_Analyses.OutdatedHpoGenesRetriever_2018-06-06.json"};		
//		args = new String[]{"E:/Groningen/Test/GeneNetwork.SubHpoPredictor/config_GeneNetwork.SubHpoPredictor_2018-03-09.json"};
		
//		args = new String[]{"E:/Groningen/Test/Analyses.ExomizerAnalysis/config_Analyses.ExomizerAnalysis_2018-07-05.json"};
//		args = new String[]{"E:/Groningen/Test/GeneNetwork.ExomizerRankRetriever/config_GeneNetwork.ExomizerRankRetriever_2018-07-09.json"};
		
		//String testClass = "GeneNetwork.CovariateAdjuster";
		//args = new String[]{"getjson",testClass,"E:/Groningen/Test/"+testClass+"/"};
			
//		String testClass = "GeneNetwork.GetGeneRanks";
		//args = new String[]{"E:/Groningen/Test/GeneNetwork.GetGeneRanks/BenchmarkSamples3/config_GeneNetwork.GetGeneRanks_2018-03-06.json"};
		//args = new String[]{"E:/Groningen/Test/GeneNetwork.GetGeneRanks/Benchmark7/config_GeneNetwork.GetGeneRanks.json"};
//		args = new String[]{"E:/Groningen/Test/GeneNetwork.UnannotatedSampleZscoreCounter/config_GeneNetwork.UnannotatedSampleZscoreCounter_2018-07-24.json"};
		
//		args = new String[]{"E:/Groningen/Test/MatrixScripts.Transpose/config_MatrixScripts.Transpose_2017-06-05.json"};
//		args = new String[]{"E:/Groningen/Test/Kallisto._Kallisto_Pipeline/config_Kallisto._Kallisto_Pipeline_2017-10-19.json"};
//		args = new String[]{"E:/Groningen/Test/GeneNetwork.WebsiteMatrixCreator/config_GeneNetwork.WebsiteMatrixCreator_HPO.json"};
		//args = new String[]{"E:/Groningen/Data/Annotation/GRCh38/PathwaysPatrick/config_MatrixScripts.GetRows_GO_F.json"};
		
		Run.run(args);
	}
}
