package MatrixScripts;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;

import Tools.FileUtils;
import Tools.Script;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.concurrent.ConcurrentCorrelation;
import umcg.genetica.math.stats.concurrent.ConcurrentCovariation;

public class RemoveDuplicates extends Script<RemoveDuplicates> 
{
	String expressionFn = "";
	String writeFn = "";
	double minCorrelationToBeDuplicate = 0.9999;
	boolean removeDuplicatesFromRows=true;//if false removes duplicates from rows
	int threads = 20;
	
	@Override
	public void run() 
	{
		try
		{
			log("Loading data");
			DoubleMatrixDataset<String, String> doubleMatrixDataset = DoubleMatrixDataset.loadDoubleData(expressionFn);

			DoubleMatrixDataset<String, String> orientedMatrix = doubleMatrixDataset;
			if(!removeDuplicatesFromRows)
				orientedMatrix=doubleMatrixDataset.viewDice();
			
			
			double[][] correlation = getCorrelation(orientedMatrix, threads) ;
			
			ArrayList<String> row_To_RowName = doubleMatrixDataset.getRowObjects();
			ArrayList<String> col_To_ColName = doubleMatrixDataset.getColObjects();
			
			ArrayList<String> samplesToRemove = new ArrayList<String>();
			
			for(int r = 0; r < correlation.length;r++)
			{
				for(int c = r+1; c < correlation[r].length;c++)
				{
					if(correlation[r][c]> minCorrelationToBeDuplicate)
					{
						if(removeDuplicatesFromRows)
						{
							String rowName = row_To_RowName.get(c);
							doubleMatrixDataset.getHashRows().remove(rowName);
							samplesToRemove.add(rowName);
						}
						else
						{
							String rowName = col_To_ColName.get(c);
							doubleMatrixDataset.getHashCols().remove(rowName);
							samplesToRemove.add(rowName);
						}
					}
				}
			}
			

			log("Writing leftover matrix to: \t" + writeFn);
			doubleMatrixDataset.save(writeFn);
			
		}catch (Exception e){e.printStackTrace();}
	}
	public double[][] getCorrelation(DoubleMatrixDataset<String, String> expression, int threads) throws IOException 
	{
		double[][] matrix = null;
		
		//use proper library to get correlations
		log("Calculating correlation");
		ConcurrentCorrelation calculator = new ConcurrentCorrelation(threads);
		matrix = calculator.pairwiseCorrelation(expression.getMatrix().toArray());
		return matrix;
	}
	
}
