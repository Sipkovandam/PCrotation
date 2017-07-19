package Tools;

import java.io.IOException;

import MatrixScripts.MatrixStruct;
import MatrixScripts.MyMatrix;

public class MatrixStats extends Script<MatrixStats>{

	public String matrixFN = "E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/TESTexpression.txt";
	public String writeFN_Rows = null;
	public String writeFN_Cols = null;
	
	public void run()
	{
		try
		{
			//rows
			MatrixStruct matrix = new MatrixStruct(matrixFN);
			MatrixStruct averagesRows = matrix.getAveragesPerRow();
			MatrixStruct stdevsRows = matrix.stDevRows();
			MatrixStruct mergeRows = averagesRows.mergeColumns(stdevsRows);
					
			if(writeFN_Rows == null)
				if(mergeRows.getRowHeaders()[0].contains("ENSG0") || mergeRows.getRowHeaders()[0].contains("ENST0"))
				{
					writeFN_Rows=FileUtils.replaceEnd(matrixFN, "_Gene_stats.txt");
					writeFN_Cols=FileUtils.replaceEnd(matrixFN, "_Sample_stats.txt");
				}
				else
				{
					writeFN_Rows=FileUtils.replaceEnd(matrixFN, "_Row_stats.txt");
					writeFN_Cols=FileUtils.replaceEnd(matrixFN, "_Col_stats.txt");
				}
	
			mergeRows.write(writeFN_Rows);
			
			//cols
			MatrixStruct colAverages = matrix.getAveragesPerCol();
			MatrixStruct colSstdevs = matrix.stDevCols();
			MatrixStruct colMerge = colAverages.mergeColumns(colSstdevs);
			colMerge.write(writeFN_Cols);
		}catch(Exception e)
		{
			e.printStackTrace();
		}
	}
}
