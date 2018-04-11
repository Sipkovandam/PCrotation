package MatrixScripts;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.lang.Math;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.stat.inference.WilcoxonSignedRankTest;

import RowAnalyses.Row;
import Tools.FileUtils;
import Tools.Script;

import org.apache.commons.math3.stat.inference.TTest;
import org.apache.commons.math3.stat.StatUtils;


public class ZscoreCalculator extends Script<ZscoreCalculator>
{
	//Converts values per row (gene) to Z-scores per gene across all samples
	public String sampleFn = null;
	public String writeFn = null;
	public String statsFn = null;
	
	@Override
	public void run()
	{
		try
		{
			MyMatrix test=new MyMatrix();
			
			
			System.out.println("Done! File written to: " + writeFn);
			if(writeFn==null)
			{
				writeFn = FileUtils.replaceEnd(sampleFn, "_zScores.txt.gz");
			}
			MyMatrix scoreStats = null;
			if(statsFn!=null)
				scoreStats=new MyMatrix(statsFn);
			
			BufferedReader sampleReader = FileUtils.createReader(sampleFn);
			BufferedWriter zscoreWriter = FileUtils.createWriter(writeFn);
			String header = sampleReader.readLine();
			zscoreWriter.write(header+"\n");
			String line = null;
			while((line=sampleReader.readLine())!=null)
			{
				Row row = Row.readRow(line);
				zscoreWriter.write(row.getRowName());
				
				if(!scoreStats.getRowHash().containsKey(row.getRowName()))
				{
					log("Error: statsFile missing rowname: " + row.getRowName() +  "\n exiting");
					System.exit(2);
				}
				
				int scoreStatIndex = scoreStats.getRowHash().get(row.getRowName());
				double avg = scoreStats.values[scoreStatIndex][0];
				double stdev = scoreStats.values[scoreStatIndex][1];;
				for(int v = 0; v < row.getValues().length; v++)
				{
					double zscore = (row.getValues()[v]-avg)/stdev;
					zscoreWriter.write("\t"+zscore);
				}
				
				zscoreWriter.write("\n");
			}
			zscoreWriter.close();
			sampleReader.close();
		}catch(Exception e){e.printStackTrace();}
		
	}

	public String getSampleFn()
	{
		return sampleFn;
	}

	public void setSampleFn(String sampleFn)
	{
		this.sampleFn = sampleFn;
	}

	public String getWriteFn()
	{
		return writeFn;
	}

	public void setWriteFn(String writeFn)
	{
		this.writeFn = writeFn;
	}

	public String getStatsFn()
	{
		return statsFn;
	}

	public void setStatsFn(String statsFn)
	{
		this.statsFn = statsFn;
	}

// oude zooi
//	public MyMatrix changeToZscores(MyMatrix samples) 
//	{
//		return changeToZscores(samples, null);
//	}
//	
//	//calculates the z-scores over the rows
//	public MyMatrix changeToZscores(MyMatrix samples, String statsFileName) 
//	{
//		System.out.println(" Memory usage:" + Runtime.getRuntime().totalMemory()/1024/1024/1024 + " GB");
//		
//		MyMatrix scoreStats = new MyMatrix(samples.rows(),2);
//		scoreStats.setRowHeaders(samples.getRowHeaders());
//		scoreStats.setColHeaders(new String[]{"Average","StDev"});
//		if(statsFileName != null)
//			scoreStats = new MyMatrix(statsFileName);
//		
//		for(int r = 0; r < samples.getRowHeaders().length; r++)//each has the values for 1 row (for all samples)
//		{
//			double mean = 0,stDev = 0;
//			if(statsFileName == null)
//			{
//				double[] rowValues = samples.getRowValues(r);
//				mean = org.apache.commons.math3.stat.StatUtils.mean(rowValues);
//				double variance = org.apache.commons.math3.stat.StatUtils.populationVariance(rowValues);
//				stDev = java.lang.Math.pow(variance,0.5);
//				scoreStats.matrix.set(r, 0, mean);
//				scoreStats.matrix.set(r, 1, stDev);
//			}
//			else
//			{
//				String geneName = samples.getRowHeaders()[r];
//				if(!scoreStats.getRowHash().containsKey(geneName))
//					System.out.println("Warning this gene is missing from the stats file:" + geneName);
//				int scoreRow = scoreStats.getRowHash().get(geneName);
//				mean = scoreStats.matrix.get(scoreRow, 0);
//				stDev = scoreStats.matrix.get(scoreRow, 1);
//			}
//			
//			for(int c = 0; c < samples.getColHeaders().length; c++)
//			{
//				double currentVal = samples.matrix.get(r, c);
//				samples.matrix.set(r, c, zScore(currentVal,mean,stDev));
//			}
//			
////			Arrays.stream(rowValues)
////				  .map(value -> zScore(value,mean,stDev))
////				  .forEach(zScore -> storePCinfo(zScore,output));
//			
//		}
//		return scoreStats;
//	}
//
//	private double zScore(double value, double mean, double stDev) {
//		if(stDev == 0) return 0;
//		return (value-mean)/stDev;
//	}
//	
//	public MyMatrix zScores(String writeFolder, String zscoreFN, MyMatrix sampleStruct, String inputAvgStdevFolder, String statsFN, boolean addStatColumns) throws IOException 
//	{
//		JuhaPCA.PCA.log("Calculating zScores");
//		MyMatrix zScoreMatrix = sampleStruct.copy();
//		MyMatrix avgStdev = null;
//		if(inputAvgStdevFolder == null)//use input stats if supplied
//			avgStdev = this.changeToZscores(zScoreMatrix, null);
//		else
//			avgStdev = this.changeToZscores(zScoreMatrix, inputAvgStdevFolder+statsFN);
//				
//		JuhaPCA.PCA.log("Writing zScores");
//		zScoreMatrix=avgStdev.mergeColumns(zScoreMatrix);
//		if(inputAvgStdevFolder !=null)
//		{
//			zScoreMatrix.setColHeader(0, "average_" + inputAvgStdevFolder);
//			zScoreMatrix.setColHeader(1, "stDev_" + inputAvgStdevFolder);
//		}
//		zScoreMatrix.write(writeFolder+zscoreFN);
//		if(inputAvgStdevFolder == null)
//			avgStdev.write(writeFolder+statsFN);
//		
//		if(addStatColumns)
//		{
//			sampleStruct = avgStdev.mergeColumns(sampleStruct);
//			if(inputAvgStdevFolder !=null)
//			{
//				sampleStruct.setColHeader(0, "average_" + inputAvgStdevFolder);
//				sampleStruct.setColHeader(1, "stDev_" + inputAvgStdevFolder);
//			}
//		}
//		return sampleStruct;
//	}
//	
}

//double zScore = 0;
//new WilcoxonSignedRankTest().wilcoxonSignedRank(arg0, arg1)
//WilcoxonSignedRankTest test = new WilcoxonSignedRankTest();
//TTest tTest = new TTest();
//tTest.
//tTest.
//test.
//samples.getRowValues(r);
