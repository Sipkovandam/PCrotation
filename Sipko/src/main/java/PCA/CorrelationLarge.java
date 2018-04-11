package PCA;

import java.io.IOException;
import java.util.ArrayList;

import MatrixScripts.MyMatrix;
import Tools.FileUtils;
import Tools.Script;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.concurrent.ConcurrentCorrelation;
import umcg.genetica.math.stats.concurrent.ConcurrentCovariation;

public class CorrelationLarge extends Script<CorrelationLarge>{
	/**
	 * 
	 */
	private static final long serialVersionUID = 5764301796123217044L;
	//creates a correlation matrix without the limitation of a java array size. Limitation is now the memory of the machine
	private String expressionFN = "E:/Groningen/Data/Test/pre_Correlation_Or_Covariance.txt.gz";
	private String writeFn = null;
	private int threads = 20;
	private boolean correlation = false;
	private boolean overRows = true;
	private transient boolean writeParralel=false;//if true writes the output in parallel proceeding to the next step with the resulting matrix

	@Override
	public void run() 
	{
		this.run(null, false);
	}
	
	public MyMatrix run(MyMatrix passMatrix, boolean writeParallel) 
	{
		MyMatrix covMat=null;
		try
		{
			writeFn = checkWriteFn();

//			int rowsWithoutVariance=expression.removeNoVariance();
//			if(rowsWithoutVariance>0)
//				FileUtils.addBeforeExtention(writeFn, "RowsWithoutVarianceRemoved");
				
			covMat = correlation(writeFn, expressionFN,correlation,threads);
			log("Writing file from which duplicates are removed:\t" + writeFn);
			covMat.write(writeFn,false,writeParallel);
			if(!writeParallel)
				System.out.println("Done, file written to: "+ writeFn);
		} catch (IOException e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
			log("Exiting");
			System.exit(2);
		}
		return covMat;
	}

	private String checkWriteFn()
	{
		if(writeFn == null)
		{
			if(correlation)
				writeFn = expressionFN.replace(".txt", "").replace(".gz", "") + "_correlation.txt.gz";
			else//if covariance
				writeFn = expressionFN.replace(".txt", "").replace(".gz", "") + "_covariance.txt.gz";
		}
		return writeFn;
	}

	public MyMatrix correlation(String writeFile, String expressionFn, boolean correlation, int threads) throws IOException 
	{
		double[][] matrix = null;
		log("Loading matrix");
		DoubleMatrixDataset<String, String> doubleMatrixDataset = DoubleMatrixDataset.loadDoubleData(expressionFn);
		
		DoubleMatrixDataset<String, String> orientedMatrix = doubleMatrixDataset;
		if(!overRows)
			orientedMatrix=doubleMatrixDataset.viewDice();
		
		
		if(correlation)
		{
			log("Calculating correlation");
			matrix = getCorrelation(orientedMatrix, threads);			
		}
		else
		{
			log("Calculating covariation");
			matrix = getCovariance(orientedMatrix, threads);	
		}
		ArrayList<String> rowNames = orientedMatrix.getRowObjects();
		MyMatrix covMat = new MyMatrix(rowNames.toArray(new String[rowNames.size()]), rowNames.toArray(new String[rowNames.size()]), matrix);
		
		return covMat;
	}
	
	public static double[][] getCorrelation(DoubleMatrixDataset<String, String> expression, int threads) throws IOException 
	{
		double[][] matrix = null;
		
		//use proper library to get correlations
		ConcurrentCorrelation calculator = new ConcurrentCorrelation(threads);
		matrix = calculator.pairwiseCorrelation(expression.getMatrix().toArray());
		return matrix;
	}
	
	public static double[][] getCovariance(DoubleMatrixDataset<String, String> expression, int threads) throws IOException 
	{
		double[][] matrix = null;
		
		//use proper library to get correlations
		ConcurrentCovariation calculator = new ConcurrentCovariation(threads);
		matrix = calculator.pairwiseCovariation(expression.getMatrix().toArray());
		return matrix;
	}

	public String getExpressionFN()
	{
		return expressionFN;
	}

	public void setExpressionFN(String expressionFN)
	{
		this.expressionFN = expressionFN;
	}

	public int getThreads()
	{
		return threads;
	}

	public void setThreads(int threads)
	{
		this.threads = threads;
	}

	public boolean isCorrelation()
	{
		return correlation;
	}

	public void setCorrelation(boolean correlation)
	{
		this.correlation = correlation;
	}

	public String getWriteFn()
	{
		return writeFn;
	}

	public void setWriteFn(String writeFn)
	{
		this.writeFn = writeFn;
	}

	public boolean isWriteParralel()
	{
		return writeParralel;
	}

	public void setWriteParralel(boolean writeParralel)
	{
		this.writeParralel = writeParralel;
	}

	public boolean isOverRows()
	{
		return overRows;
	}

	public void setOverRows(boolean overRows)
	{
		this.overRows = overRows;
	}
}