package Tests;

import Tools.Run;

public class RunTest
{
	
	public static void main(String[] args) throws Exception
	{
		//args = new String[]{"getjson","Tests.Test","E:/Groningen/Test/Run/"};
		//args = new String[]{"E:/Groningen/Test/Run/config_Tests.Test_2017-01-18.json"};
		//args = new String[]{"E:/Groningen/Test/JSON/STAR._STAR_Pipeline_config_2017-01-16_initiated.json"};
		//args = new String[]{"getjson","STAR._STAR_Pipeline","E:/Groningen/Splicing/"};
		//args = new String[] { "getjson", "PCA.PCApipeline", "E:/Groningen/Test/JSON/" };
		//args = new String[] { "E:/Groningen/Test/JSON/PCApipeline/config_PCA.PCApipeline_2017-02-10.json" };
		
		String testClass = "R.SplicingAnalysisPipeline";
		//args = new String[]{"getjson",testClass,"E:/Groningen/Test/"+testClass+"/"};
		args = new String[]{"E:/Groningen/Test/"+testClass+"/config_R.SplicingAnalysisPipeline_2017-06-25.json"};

		
		//args = new String[]{"getjson", "Tests.Test", "E:/Groningen/Test/Tests.Test/"};
		//args = new String[]{"E:/Groningen/Test/Tests.Test/config_Tests.Test_2017-06-21.json"};
		
		Run.run(args);
	}
}
