package MatrixScripts;

import Tools.Script;

public class NoVarianceRowRemover extends Script<NoVarianceRowRemover>
{
	String fn = null;
	String writeName = null;
	
	@Override
	public void run()
	{
		MyMatrix matrix = new MyMatrix(fn);
		matrix.removeNoVariance();
	}	
}
