package PCA;

import java.io.File;
import java.io.IOException;

import MatrixScripts.MyMatrix;
import Tools.FileUtils;
import Tools.Script;
import eqtlmappingpipeline.normalization.NormalizationConsoleGUI;

public class QuantileNormalize extends Script<QuantileNormalize>
{
	//Quantile normalizes a dataset
	String expressionFn = "";
	String writeFolder = "";
	
	@Override
	public void run() 
	{
		try
		{
			FileUtils.makeDir(writeFolder);
			new NormalizationConsoleGUI(new String[]{"--qqnorm", "--in", expressionFn, "--out", writeFolder});
		}catch(Exception e){e.printStackTrace();}
	}

	public void quantileNormalize(	String expressionStruct,
											String writeFolder)
	{
		this.expressionFn=expressionStruct;
		this.writeFolder=writeFolder;
		this.run();
	}

//	public void quantileNormalize(MyMatrix expressionStruct) throws IOException {		
//		
//	}
}
