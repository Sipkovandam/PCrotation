package STAR;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import MatrixScripts.RowStatisticsGetter;
import MatrixScripts.CenterPerRow;
import MatrixScripts.LaneMerger;
import MatrixScripts.MyMatrix;
import MatrixScripts.Transpose;
import PCA.CorrelationLarge;
import PCA.Pca;
import PCA.PcaPipelineLite;
import Slurm.Slurm;
import TextEditing.FileSplitter;
import Tools.FileUtils;
import Tools.Runnable;
import Tools.Script;

public class ExonCorrectionPipeline extends Script<ExonCorrectionPipeline>
{
	private String exonExpressionMergedFn = null;
	private String writeFolder = null;
	private int nThreads = 8;
	private boolean spliceInsteadOfExon = false;
	
	@Override
	public void run()
	{
		try
		{
			List<Runnable> steps = initiate();

			// run the selected steps
			MyMatrix passMatrix = null;
			for (Runnable step : steps)
				passMatrix=step.run(passMatrix);
		} catch (Exception e)
		{
			e.printStackTrace();
		}
	}
	
	private List<Runnable> initiate() throws CloneNotSupportedException
	{
		FileUtils.makeDir(writeFolder);
					
		//FileSplitter;
		String writeFolderUncorrectedPerGene = writeFolder + "/expressionUncorrectedPerGene/";
		FileSplitter fileSplitter = new FileSplitter();
		fileSplitter.setFn(exonExpressionMergedFn);
		fileSplitter.setGzipFiles(false);
		fileSplitter.setSpliceNamesCol(3);//not relevant actually I think
		fileSplitter.setSplitRegex("__");
		fileSplitter.setSplitFirstColumnInstead(true);
		
		if(spliceInsteadOfExon)
		{
			fileSplitter.setSpliceNamesCol(3);
			fileSplitter.setSplitRegex(",");
			fileSplitter.setSplitFirstColumnInstead(false);
		}
		fileSplitter.setWriteFolder(writeFolderUncorrectedPerGene);

		//RatioCalculator; calculate relative usage of each splice site for each gene
		RatioCalculator ratioCalculator = new RatioCalculator();
		ratioCalculator.setExpressionFolder(fileSplitter.getWriteFolder());
		ratioCalculator.setWriteFn(new File(ratioCalculator.getExpressionFolder()).getParent()+"/ratios.txt.gz");
		
		
		//Transpose; transpose the matrix so the slice variants are on the columns
		Transpose transpose = new Transpose();
		transpose.setFileName(ratioCalculator.getWriteFn());
		transpose.setWriteFn(FileUtils.removeExtention(transpose.getFileName())+ "_transposed.txt.gz");
		
		String pcaDir = FileUtils.removeExtention(transpose.getFileName())+"Pca/";
		new File(pcaDir).mkdirs();
		
		//removeRows withou
		
		//CenterRows
		CenterPerRow centerPerRow = new CenterPerRow();
		centerPerRow.setFileName(transpose.getWriteFn());
		centerPerRow.setAveragesWriteFn(FileUtils.removeExtention(transpose.getWriteFn())+"_sampleAverages.txt");
		centerPerRow.setWriteFn(pcaDir+"MATRIX_Centered.txt.gz");
		
		//CorrelationLarge; calculate a covariance matrix over the samples dimension in the ratios file;
		CorrelationLarge correlationLarge = new CorrelationLarge();
		correlationLarge.setExpressionFN(centerPerRow.getWriteFn());
		correlationLarge.setCorrelation(false);
		correlationLarge.setThreads(nThreads);
		correlationLarge.setWriteFn(pcaDir+"ratios_rowCentered_covariance.txt.gz");
		
		//run PCA over the covariance matrix
		Pca pca = new Pca();//this requires the whole matrix to be in 1 array.
		pca.setInputMatrix(correlationLarge.getWriteFn());
		pca.setEigenVectorWriteFn(pcaDir+ "SAMPLE.eigenvectors.txt");
		pca.setEigenValueWriteFn(pcaDir+"SAMPLE.eigenvalues.txt");
		
		//calculate the averages per gene (not used for calculations, just to indicate which genes should be included)
		RowStatisticsGetter averagesPerRow = new RowStatisticsGetter();
		averagesPerRow.setFileName(ratioCalculator.getWriteFn());
		averagesPerRow.setWriteFn(pcaDir+"SAMPLE_Norm_GeneAverages.txt");
		averagesPerRow.setAbsolute(false);
		
		//run PC correction on the ratios file
		String correctedDir = pcaDir+"ratiosCorrected/";
		PcaPipelineLite PCApipelineLite = new PcaPipelineLite();
		PCApipelineLite.setSampleFile(ratioCalculator.getWriteFn());
		PCApipelineLite.setRunMode(2);
		PCApipelineLite.setPcaOverGenes(false);//not used, it uses centered matrix from PCA
		PCApipelineLite.setCenterGenes(false);//not used, it uses centered matrix from PCA
		PCApipelineLite.setCenterSamples(true);//not used, it uses centered matrix from PCA
		PCApipelineLite.setCorrectResultsForSTdevs(false);//not used, it uses centered matrix from PCA
		PCApipelineLite.setWriteFolder(pcaDir);
		PCApipelineLite.setLog2(false);//not used, it uses centered matrix from PCA
		PCApipelineLite.setDeSeqNorm(false);//not used, it uses centered matrix from PCA
		PCApipelineLite.setWriteFolderCorrected(correctedDir);
		PCApipelineLite.setPCs("1-100");
		
		FileSplitter fileSplitterCorrectedRatios = new FileSplitter();
		fileSplitterCorrectedRatios.setFn(PCApipelineLite.getWriteFolderCorrected()+"PC_"+PCApipelineLite.getPCs()+".txt.gz");
		fileSplitterCorrectedRatios.setGzipFiles(false);
		fileSplitterCorrectedRatios.setSplitRegex("__");
		fileSplitterCorrectedRatios.setSplitFirstColumnInstead(true);
		fileSplitterCorrectedRatios.setWriteFolder(correctedDir+"perGene/");

		// add steps (in right order of course ;))
		List<Runnable> steps = new ArrayList<>();
		steps.add(fileSplitter);
		steps.add(ratioCalculator);
		steps.add(transpose);
		steps.add(centerPerRow);
		steps.add(correlationLarge);
		steps.add(pca);
		steps.add(averagesPerRow);
		steps.add(PCApipelineLite);
		steps.add(fileSplitterCorrectedRatios);
		return steps;

	}
}
