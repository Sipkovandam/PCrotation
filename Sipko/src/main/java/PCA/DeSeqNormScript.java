package PCA;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import MatrixScripts.MyMatrix;
import Tools.FileUtils;
import Tools.Script;

public class DeSeqNormScript extends Script<DeSeqNormScript>
{
	//This is doing the DESeq normalization, which does not use replicate information.
	
	String expressionFN = "E:/Groningen/Scripts/Tests/Rlog.java/Samples.txt";
	
	String writeFolder = null;//if null becomes new File(expressionFN).getParent()+"/";
	String writeFn = null;//if null becomes new File(expressionFN).getParent()+"/";
	String geoFN = null;//if null calculates geometric means based on this dataset
	boolean writeAll = true;	//write all intemediary files too
	double logAdd = 1;//previously used to multiply the results by this number, but seems pointless since it does not have any effect, so does not do anything anymore
	boolean log = false;
	boolean roundValues = false;//rounds expression values to whole counts
	
	@Override	
	public void run() 
	{
		try
		{
			init();
			
			MyMatrix expression = new MyMatrix(expressionFN);
			
			if(roundValues)
				expression.roundValues();
			double start = System.nanoTime();
			rLog(writeFolder, expression, writeAll, geoFN);
			double end = System.nanoTime();
			System.out.println((end-start)/1000/1000/1000 + " sec");
			if(log)
			{
				expression.log2Transform(logAdd);
				//write the results
					writeFn= FileUtils.addBeforeExtention(writeFn, ".Log2");
				expression.write(this.writeFn);
			}
		}catch(Exception e){e.printStackTrace();}
	}

	private void init()
	{
		DeSeqNorm.expressionFN=this.expressionFN;
		DeSeqNorm.writeFolder=this.writeFolder;
		DeSeqNorm.geoFN=this.geoFN;
		DeSeqNorm.writeAll=this.writeAll;
		DeSeqNorm.logAdd=this.logAdd;
		DeSeqNorm.log=this.log;
		DeSeqNorm.roundValues=this.roundValues;	
		if(writeFolder == null)
			writeFolder = new File(expressionFN).getParent()+"/DESeqNorm/";
		

		if(writeFn==null)
			writeFn =  writeFolder + new File(expressionFN).getName().replace(".txt", "").replace(".gz","") + ".DESeqNorm.txt.gz";
		
		if(!new File(writeFolder).exists())
			new File(writeFolder).mkdirs();
		
	}
	

	public void rLog(String writeFolder, MyMatrix expressionStruct, boolean writeAll, String geoFN) throws IOException 
	{
		String swapFN = writeFolder + "swapFile.txt";
		expressionStruct.write(swapFN);
		JuhaPCA.PCA.log(" 6. Deseq normalization");
		if(geoFN!= null && new File(geoFN).exists())
			rLog(expressionStruct, writeFolder, swapFN, new MyMatrix(geoFN), null);
		else
			rLog(expressionStruct, writeFolder, swapFN, null, geoFN);
		if(writeAll)
			expressionStruct.write(this.writeFn);
	}
	
	public void rLog(MyMatrix expressionStruct, String writeFolder, String fileName, String writeGeoFN) throws IOException 
	{
		rLog(expressionStruct, writeFolder, fileName, null, writeGeoFN);
	}
	
	private MyMatrix getGeoMeans(MyMatrix expressionStruct, String writeFolder, String writeGeoFN) throws IOException {
		//log transform data so values do not become to large (and you save a lot of speed (you can add instead of multiply))
		expressionStruct.logTransform(10,0);
		//create the matrix where the results are stored
		MyMatrix geoMean = new MyMatrix(expressionStruct.rows(),1);
		geoMean.setRowHeaders(expressionStruct.getRowHeaders());
		//Creates a matrix where you can put the genes you want to keep (would normally use a hashtable, but this works with MatrixStruct.keeprows called later)
		MyMatrix keepGenes = new MyMatrix(expressionStruct.rows(),1);
		for(int r =0; r < expressionStruct.rows(); r++)
		{
			//Calculate the geometric mean (based on the logged values, which is why it is summing the row instead of multiplying. Also why it is dividing instead of ^(1/x))
			double gM = expressionStruct.sumRow(r)/expressionStruct.cols();
			if(Double.isFinite(gM))
			{
				geoMean.matrix.set(r,0,gM);
				geoMean.matrix.set(r,0,Math.pow(10,geoMean.matrix.get(r,0)));//raise 10^X to unlog the resulting values again
				keepGenes.setRowHeader(r, expressionStruct.getRowHeaders()[r]);
			}
		}
		
		geoMean.keepRows(keepGenes); 
		if(writeGeoFN==null)
			writeGeoFN = writeFolder+"geoMean.txt";
		System.out.println("geofilename=" + writeGeoFN);
		geoMean.write(writeGeoFN);
		
		return geoMean;
	}
	
	public void rLog(MyMatrix expression, String writeFolder, String fileName, MyMatrix geoMean, String writeGeoFN) throws IOException 
	{
		if(geoMean == null)
		{
			geoMean = getGeoMeans(expression, writeFolder, writeGeoFN);
			//need to read the matrix again here because we logged it before...
			expression.readFile(fileName);
		}
		MyMatrix denominators = getDenominators(expression, fileName, geoMean, writeFolder);
		
		expression.readFile(fileName);//read the file again since we kept only the rows that have no 0 values.
		
		//calculate the normalized readcounts
		for(int c = 0; c < expression.cols(); c++)
			for(int r = 0; r < expression.rows(); r++)
				expression.matrix.set(r,c,expression.matrix.get(r,c)/denominators.matrix.get(c,0));	
	}
	
	private MyMatrix getDenominators(MyMatrix expressionStruct, String fileName, MyMatrix geoMean, String writeFolder) throws IOException {
		geoMean.keepRows(expressionStruct);//keep only the rows that are also in the geomean file
		expressionStruct.write(fileName.replace(".txt", "geoRowsOnly.txt"));
		MyMatrix denominators = new MyMatrix(expressionStruct.cols(),1);
		denominators.setRowHeaders(expressionStruct.getColHeaders());
		
		//determine the denominator per sample
		for(int c = 0; c < expressionStruct.cols(); c++)
		{
			//get the correct denominator
			double[] column = new double[expressionStruct.rows()];
			for(int r = 0; r < column.length; r++)
			{
				column[r] = (expressionStruct.matrix.get(r,c)) / geoMean.matrix.get(r,0);
			}
			Arrays.sort(column);
			
			org.apache.commons.math3.stat.descriptive.rank.Median med = new org.apache.commons.math3.stat.descriptive.rank.Median();
			denominators.matrix.set(c,0,med.evaluate(column));
		}
		if(writeFolder != null)
			denominators.write(writeFolder + "Denominators.txt");
		return denominators;
	}

	public String getExpressionFN()
	{
		return expressionFN;
	}

	public void setExpressionFN(String expressionFN)
	{
		this.expressionFN = expressionFN;
	}

	public String getWriteFolder()
	{
		return FileUtils.makeFolderNameEndWithSlash(writeFolder);
	}

	public void setWriteFolder(String writeFolder)
	{
		this.writeFolder = writeFolder;
	}

	public String getGeoFN()
	{
		return geoFN;
	}

	public void setGeoFN(String geoFN)
	{
		this.geoFN = geoFN;
	}

	public double getLogAdd()
	{
		return logAdd;
	}

	public void setLogAdd(double logAdd)
	{
		this.logAdd = logAdd;
	}

	public boolean isLog()
	{
		return log;
	}

	public void setLog(boolean log)
	{
		this.log = log;
	}

	public boolean isRoundValues()
	{
		return roundValues;
	}

	public void setRoundValues(boolean roundValues)
	{
		this.roundValues = roundValues;
	}

	public boolean isWriteAll()
	{
		return writeAll;
	}

	public void setWriteAll(boolean writeAll)
	{
		this.writeAll = writeAll;
	}

	public String getWriteFn()
	{
		return writeFn;
	}

	public void setWriteFn(String writeFn)
	{
		this.writeFn = writeFn;
	}
}
