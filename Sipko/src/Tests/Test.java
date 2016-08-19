package Tests;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import org.junit.rules.TemporaryFolder;

import PCA.Matrix;
import PCA.MatrixStruct;
import PCA.RLog;
import PCA_testcases.CompareFiles;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

public class Test {
	
	public static void main(String[] args) throws IOException
	{	
//		String geoFN= "E:/Groningen/Scripts/Tests/Rlog.java/DESeqNorm/Samples.DESeqNorm.txt.gz";
//		CompareFiles.compare(geoFN,"TestCaseFiles/Result_Samples.txt",true);
//		
		String fileName = "TestCaseFiles/Samples.txt";
		String folderName = new File(fileName).getParent()+"/DESeqNorm/";
//		
//		TemporaryFolder folder = new TemporaryFolder();
//		System.out.println(folderName);
//		File tempFolder = folder.newFolder("test");
		
		TemporaryFolder testFolder = new TemporaryFolder();
		 File tempFile = testFolder.newFile("file.txt");
	        File tempFolder = testFolder.newFolder("folder");
	        System.out.println("Test folder: " + testFolder.getRoot());
	        // test...
		
		//		folder.newFolder(folderName);
//		File f = new File(folderName);
//		System.out.println(f.exists());
//		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(folderName+"test.txt")));
//		writer.write("test");
//		writer.close();
	}

}
