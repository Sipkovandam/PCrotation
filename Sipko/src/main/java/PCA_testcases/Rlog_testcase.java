package PCA_testcases;

import java.io.IOException;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import MatrixScripts.MatrixStruct;
import MatrixScripts.MyMatrix;
import PCA.DeSeqNorm;

public class Rlog_testcase {

	@Rule
	public TemporaryFolder folder = new TemporaryFolder();
	@Test
	public void test() throws IOException 
	{

		System.out.println("Test folder: " + folder.getRoot());
		
		String fileName = "TestCaseFiles/Samples.txt";
		String writeFolder = folder.getRoot().getAbsolutePath()+"/";
		System.out.println(folder.getRoot().exists());
		
		MyMatrix testFile = new MyMatrix(fileName);
		DeSeqNorm.rLog(testFile, writeFolder, fileName, null);
		String rlogFN= writeFolder+"Samples.DESeqNorm.txt.gz";
		testFile.write(rlogFN);
		
		String geoFN= writeFolder+"/geoMean.txt";
		CompareFiles.compare(geoFN,"TestCaseFiles/Result_Geomean.txt",true);
		String denominatorsFN= writeFolder+"Denominators.txt";	
		CompareFiles.compare(denominatorsFN,"TestCaseFiles/Result_Denominators.txt",true);
		CompareFiles.compare(rlogFN,"TestCaseFiles/Result_Samples.txt",true);
	}

}
