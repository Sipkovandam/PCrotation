package GeneNetwork;

import MatrixScripts.GetCols;
import MatrixScripts.MyMatrix;
import PCA.CorrelationLarge;
import PCA.DeSeqNormScript;
import Tools.FileUtils;
import Tools.Script;

public class GeneNetwork_Pipeline extends Script<GeneNetwork_Pipeline>
{
	String rawDataFn = "";
	String samplesToInclude = "";
	String covariateMatrixFn = "";
	
	int nThreads = 20;
	
	@Override
	public void run()
	{
		GetCols colGetter = new GetCols();
		colGetter.setFileName(this.rawDataFn);
		colGetter.setFileName2(samplesToInclude);
		colGetter.run();
		
		DeSeqNormScript deSeqNorm = new DeSeqNormScript();
		deSeqNorm.setExpressionFN(colGetter.getWriteName());
		deSeqNorm.setWriteAll(false);
		deSeqNorm.setLog(true);
		deSeqNorm.setLogAdd(1);
		deSeqNorm.run();
		
		CovariateAdjuster covariateAdjuster = new CovariateAdjuster();
		covariateAdjuster.setExpressionMatrix(deSeqNorm.getWriteFn());
		covariateAdjuster.setCovariateMatrix(covariateMatrixFn);
		
		String expressionFN = FileUtils.addBeforeExtention(covariateMatrixFn,".ProbesWithZeroVarianceRemoved.CovariatesRemoved");
		CorrelationLarge correlationLarge = new CorrelationLarge();
		correlationLarge.setExpressionFN(expressionFN);
		correlationLarge.setOverRows(true);
		correlationLarge.setThreads(nThreads);
		correlationLarge.setWriteFn(FileUtils.addBeforeExtention(expressionFN,".Correlation"));
		
	}
}
