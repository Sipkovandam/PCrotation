package PCA;

import java.io.IOException;

import MatrixScripts.MyMatrix;
import Tools.Script;
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
	private transient boolean writeParralel=false;//if true writes the output in parallel proceeding to the next step with the resulting matrix
	
	public void run() 
	{
		this.run(null, false);
	}
	
	public MyMatrix run(MyMatrix passMatrix, boolean writeParallel) 
	{
		//disabled this, because the next function has to accept this matrix or it won't work and you never know which script comes after this. Need to find solution
		MyMatrix covMat=null;
		try
		{
			MyMatrix expression = passMatrix;
			if(passMatrix==null)
				expression = new MyMatrix(expressionFN);

			writeFn = checkWriteFn();

			covMat = correlation(writeFn, expression,correlation,threads);
			covMat.write(writeFn,false,writeParallel);
			if(!writeParallel)
				System.out.println("Done, file written to: "+ writeFn);
		} catch (IOException e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
			p("Exiting");
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

	public static MyMatrix correlation(String writeFile, MyMatrix expression, boolean correlation, int threads) throws IOException 
	{
		double[][] matrix = null;
		if(correlation)
		{
			ConcurrentCorrelation calculator = new ConcurrentCorrelation(threads);
			matrix = calculator.pairwiseCorrelation(expression.values);			
		}
		else
		{
			ConcurrentCovariation calculatorGenes = new ConcurrentCovariation(threads);
			matrix = calculatorGenes.pairwiseCovariation(expression.values);
		}
		MyMatrix covMat = new MyMatrix(expression.rowNames, expression.rowNames, matrix);
		MyMatrix averages = covMat.getAverageCols(true);
		averages.write(writeFile.replace(".txt", "_absoluteAverages.txt"));
		
		return covMat;
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
}
