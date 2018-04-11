package GeneNetwork;

import java.io.File;

import Tools.FileUtils;
import Tools.Script;
import eqtlmappingpipeline.normalization.NormalizationConsoleGUI;

public class CovariateAdjuster extends Script<CovariateAdjuster>
{
	//takes covariates out of data by doing a PCA over teh covariate matrix and regressing out these signals
	String expressionMatrix = "";
	String covariateMatrix = "";
	String writeFolder = null;
	String sampleIncludeList = "";
	
	@Override
	public void run()
	{
		try
		{
			if(writeFolder==null)
				writeFolder=new File(expressionMatrix).getParent();
			FileUtils.makeDir(writeFolder);
			new NormalizationConsoleGUI(new String[]{"--adjustcovariates", "--cov", covariateMatrix, "--covpca", "--in", expressionMatrix, "--out", writeFolder, "--sampleInclude", sampleIncludeList});
		}catch(Exception e){e.printStackTrace();}
	}

	public String getExpressionMatrix()
	{
		return expressionMatrix;
	}

	public void setExpressionMatrix(String expressionMatrix)
	{
		this.expressionMatrix = expressionMatrix;
	}

	public String getCovariateMatrix()
	{
		return covariateMatrix;
	}

	public void setCovariateMatrix(String covariateMatrix)
	{
		this.covariateMatrix = covariateMatrix;
	}

	public String getWriteFolder()
	{
		return writeFolder;
	}

	public void setWriteFolder(String writeFolder)
	{
		this.writeFolder = writeFolder;
	}

}
