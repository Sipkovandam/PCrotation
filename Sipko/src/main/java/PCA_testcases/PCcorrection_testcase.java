package PCA_testcases;

import static org.junit.Assert.*;

import java.io.File;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import Tools.Toolkit;

public class PCcorrection_testcase {

	@Rule
	public TemporaryFolder folder = new TemporaryFolder();
	
	@Test
	public void test() throws Exception {
		String jsonFile = "TestCaseFiles/PCcorrection/eigenvectors/config.json";
		String[] args = new String[]{"PCcorrection", "json="+jsonFile};
		Toolkit.main(args);
		
		File correctFolder = new File("TestCaseFiles/PCcorrection/correctResults/");
		File[] correctFiles = correctFolder.listFiles();
		for(File correctFile : correctFiles)
		{
			String correct = correctFile.getAbsolutePath();
			String testResult = "TestCaseFiles/PCcorrection/eigenvectors/TESTexpression/"+correctFile.getName();
			System.out.println("Comparing: " + testResult + "\n" +
							   "To         "+ correct);
			CompareFiles.compare(testResult, correct);
		}
		System.out.println("Test completed");
		//	fail("Not yet implemented");
	}

}
