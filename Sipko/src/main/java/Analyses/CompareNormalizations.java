package Analyses;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import PCA.MatrixStruct;
import Tools.FileSearcher;

public class CompareNormalizations {

	public static void main(String[] args) throws IOException {
		String folderName = "E:/Groningen/Data/PublicSamples/Test13/";
		String writeFN = folderName + "significanceTests.txt";
		String checkFilesFN = folderName + "checkFiles.txt";
		String chr21FN = "E:/Groningen/Data/GenePositionInfo_Chr21Only.txt";


		FileSearcher fileSearcher = new FileSearcher();
		fileSearcher.setFolders(folderName);
		fileSearcher.setWriteName(checkFilesFN);
		fileSearcher.setSearchStrings(new String[] { "PC_1" });
		fileSearcher.setForbiddenStrings( new String[] { "Smoothed101" });
		fileSearcher.run();

		BufferedReader reader = new BufferedReader(new FileReader(new File(checkFilesFN)));
		String line = null;
		int nLines = 0;
		while ((line = reader.readLine()) != null)
			nLines++;

		MatrixStruct chr21 = new MatrixStruct(chr21FN);
		MatrixStruct results = new MatrixStruct(nLines, 2);
		results.setColHeaders(new String[] { "significance", "AUC" });

		reader = new BufferedReader(new FileReader(new File(checkFilesFN)));
		line = null;
		int r = 0;
		while ((line = reader.readLine()) != null)// each line is a file in
													// which we want to compare
													// the chr21 genes to the
													// other genes
		{
			try {
				System.out.println("line = " + line + " file= " + r + "/" + results.rows());
				if (!new File(line).exists() || line.contains(".xls")) {
					results.setRow(r, line.replace("E:\\Groningen\\Data\\PublicSamples\\Test13\\", "")
							.replace("18DownSyndrome26Normal2Cancer_counts\\", ""), new double[] { 0, 0 });
					r++;
					continue;
				}

				MatrixStruct samples = new MatrixStruct(line, -1, 19);// load
																		// first
																		// 19
																		// columns
																		// (first
																		// 19 ==
																		// down
																		// syndrome)

				if (true) {
					chr21 = new MatrixStruct(chr21FN);
					chr21.keepRows1Matrix(samples);// only keep the genes on
													// chr21 that are also in
													// the samples & order genes
													// same way
					System.out.println("sample.rows = " + samples.rows() + " chr21.rows " + chr21.rows());
				}
				MatrixStruct sample = samples.getAveragesPerRow();// get
																	// averages

				double[] onChr21 = new double[chr21.rows()];
				double[] others = new double[sample.rows() - chr21.rows()];
				int chr21Index = 0, othersIndex = 0;

				for (int gene = 0; gene < sample.rows(); gene++) {
					if (chr21.rowHash.containsKey(sample.getRowHeaders()[gene])) {
						onChr21[chr21Index] = sample.matrix.get(gene, 0);
						chr21Index++;
					} else {
						// if(line.contains("directPCA_Rlog_0.2") &&
						// line.contains("PC_1-0_DevidedBySTdevs.txt"))
						// System.out.println("sample.rows = " + sample.rows()+
						// " chr21.rows " + chr21.rows() + "gene =" + gene + "
						// sample.matrixRows= " + sample.matrix.numRows()
						// + " sample.matrixCols= " + sample.matrix.numColumns()
						// + " othersIndex " + othersIndex);
						others[othersIndex] = sample.matrix.get(gene, 0);
						othersIndex++;
					}
				}
				WilcoxonMannWhitney wmw = new WilcoxonMannWhitney();
				double pValue = wmw.returnWilcoxonMannWhitneyPValue(onChr21, others);
				double auc = wmw.getAUC();
				System.out.println("pValue =" + pValue + " AUC " + auc);
				results.setRow(r, line.replace("E:\\Groningen\\Data\\PublicSamples\\Test13\\", "")
						.replace("18DownSyndrome26Normal2Cancer_counts\\", ""), new double[] { pValue, auc });
				r++;
			} catch (Exception exception) {
				exception.printStackTrace();
				results.setRow(r, line.replace("E:\\Groningen\\Data\\PublicSamples\\Test13\\", "")
						.replace("18DownSyndrome26Normal2Cancer_counts\\", ""), new double[] { 0, 0 });
				r++;
				continue;
			}

		}
		results.write(writeFN);
		System.out.println("Done, file writen to:" + writeFN);
	}

}
