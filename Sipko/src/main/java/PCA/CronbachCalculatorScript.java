package PCA;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

import JuhaPCA.CronbachCalculator;
import Tools.Script;

public class CronbachCalculatorScript extends Script<CronbachCalculatorScript>
{
	String expressionFn=null;
	String eigenVectorFn=null;
	String scoreFn=null;

	@Override
	public void run()
	{
		try
		{
			CronbachCalculator.cronbach(Paths.get(expressionFn), Paths.get(eigenVectorFn), Paths.get(scoreFn), false);
		} catch (IOException e)
		{
			e.printStackTrace();
		}
	}
}
