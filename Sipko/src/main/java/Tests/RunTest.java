package Tests;

import Tools.Run;

public class RunTest
{
	
	public static void main(String[] args) throws Exception
	{
		//args = new String[] { "E:/Groningen/Test/JSON/PCApipeline/config_PCA.PCApipeline_2017-02-10.json" };
		
		String testClass = "Gdio.QueryConverter";
		//args = new String[]{"getjson",testClass,"E:/Groningen/Test/"+testClass+"/"};
		args = new String[]{"E:/Groningen/Test/Gdio.QueryConverter/config_Gdio.QueryConverter_2018-02-01.json"};
		
		
		//String testClass = "GeneNetwork.CovariateAdjuster";
		//args = new String[]{"getjson",testClass,"E:/Groningen/Test/"+testClass+"/"};
//		args = new String[]{"E:/Groningen/Data/GeneNetwork/Recount2/2_CovariateRemoval/config_MatrixScripts.Transpose_2017-06-05.json"};
		
		
//		args = new String[]{"E:/Groningen/Test/Gdio.MolgenisEmxParser/config_Gdio.MolgenisEmxParser_2018-01-24.json"};
		Run.run(args);
	}
}
