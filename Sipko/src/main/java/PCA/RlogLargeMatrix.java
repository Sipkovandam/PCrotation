package PCA;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import Tools.JSONutil;
import Tools.Script;

public class RlogLargeMatrix extends Script<RlogLargeMatrix>
{
	//This is doing the DESeq normalization, which does not use replicate information.
	//Also works for matrixes of larger then MAXint sizes
	
	//static String expressionFN = "E:/Groningen/Data/Annique/LLD_and_BIOS_Kallisto_EstimatedTranscriptCounts.ProbesWithZeroVarianceRemoved_removedBadSamples.txt.gz";
	String expressionFN = null;//"E:/Groningen/Data/RNAseq_clinic/5GPM/GRCh38/counts_GENES_ScaffoldGenesRemovedcounts_DiscardedRemoved.txt.gz";//"/groups/umcg-wijmenga/tmp04/umcg-jkarjalainen/31.txt";
	String writeFolder = null;//if null becomes new File(expressionFN).getParent()+"/";
	String geoFN = null;//"E:/Groningen/Scripts/Tests/Rlog.java/DESeqNorm2/geoMean.txt";//"E:/Groningen/Data/Juha/Genes31995/DuplicatesRemoved/geoMean.txt";//if null calculates geometric means based on this dataset
	boolean writeAll = true;	//write all intemediary files too
	double logAdd = 0.5;//previously used to multiply the results by this number, but seems pointless since it does not have any effect, so does not do anything anymore
	boolean roundValues = true;//rounds expression values to whole counts
	
	public RlogLargeMatrix()
	{
		
	}
	
	RlogLargeMatrix(String[] args)
	{
		checkArgs(args);
			run();
	}
	
	public void run()
	{	try 
	{	
			p("Input filename =" + expressionFN);
			if(writeFolder == null)
				writeFolder = new File(expressionFN).getParent()+"/DESeqNorm/";
			
			if(!new File(writeFolder).exists())
				new File(writeFolder).mkdir();
			
			Matrix expression = new Matrix(expressionFN);
			if(roundValues)
				expression.roundValues();
				rLog(writeFolder, expression, writeAll, geoFN);
			p(" 7. Log transforming");
			expression.logTransform(2,logAdd);//adds +0.5 before log
			expression.write(writeFolder + new File(expressionFN).getName().replace(".txt", "").replace(".gz","")+".DESeqNorm.Log2.txt.gz");
			end();
		}catch(Exception e){e.printStackTrace();}
	}
	
	private void writeVars()
	{
		//JSONutil<RlogLargeMatrix> writer = new JSONutil<RlogLargeMatrix>();
		//writer.
		writeConfig(jsonFN, this);
	}
	public String getWritePath(String name)
	{
		if(!name.contains("\\") && !name.contains("/"))
			return getFolderName(this.writeFolder)+name;
		return name;
	}
	public String getFolderName(String fn) 
	{
		if(!fn.endsWith("/") && !fn.endsWith("\\"))
			fn = fn+"/";
		return fn;
	}
	
	public void rLog(String writeFolder, Matrix expression, boolean writeAll, String geoFN) throws IOException 
	{
		MatrixStruct geoMean = null;
		if(geoFN!= null)
			geoMean = new MatrixStruct(geoFN);
		expression.putGenesOnRows();
		String swapFN = writeFolder + "swapFile.txt";
		expression.write(swapFN);
		p(" 6. Rlog without log");
		String correctedNotLogged =  writeFolder + new File(expressionFN).getName().replace(".txt", "").replace(".gz","") + ".DESeqNorm.txt.gz";
		rLog(expression, writeFolder, swapFN, geoMean, null);
		
		if(writeAll)
			expression.write(correctedNotLogged);
	}
	
	public void rLog(Matrix expression, String writeFolder, String fileName, String writeGeoFN) throws IOException 
	{
		rLog(expression, writeFolder, fileName, null, writeGeoFN);
	}
	public void rLog(Matrix expression, String writeFolder, String fileName, MatrixStruct geoMean, String writeGeoFN) throws IOException 
	{
		if(geoMean == null)
		{
			geoMean = getGeoMeans(expression, writeFolder, writeGeoFN);
			//need to read the matrix again here...
			expression.readFile(fileName);
		}

		Matrix geoMeanMat = new Matrix(geoMean);
		geoMeanMat.keepRows(expression);//keep only the rows that are also in the geomean file
		
		expression.write(fileName.replace(".txt", "geoRowsOnly.txt"));
		MatrixStruct denominators = new MatrixStruct(expression.cols(),1);
		denominators.setRowHeaders(expression.getColHeaders());
		
		//determine the denominator per sample
		for(int c = 0; c < expression.cols(); c++)
		{
			//get the correct denominator
			double[] column = new double[expression.rows()];
			for(int r = 0; r < column.length; r++)
			{
				column[r] = (expression.matrix.get(r,c)) / geoMean.matrix.get(r,0);
			}
			Arrays.sort(column);
			
			org.apache.commons.math3.stat.descriptive.rank.Median med = new org.apache.commons.math3.stat.descriptive.rank.Median();
			denominators.matrix.set(c,0,med.evaluate(column));
		}
		if(writeFolder != null)
			denominators.write(writeFolder + "Denominators.txt");
		
		expression.readFile(fileName);
		//calculate the normalized readcounts
		for(int c = 0; c < expression.cols(); c++)
		{
			for(int r = 0; r < expression.rows(); r++)
			{
				expression.matrix.set(r,c,expression.matrix.get(r,c)/denominators.matrix.get(c,0));	
			}
		}
	}
	
	private MatrixStruct getGeoMeans(Matrix expression ,String writeFolder, String writeGeoFN) throws IOException {
		expression.logTransform(10,0);
		MatrixStruct geoMean = new MatrixStruct(expression.rows(),1);
		geoMean.setRowHeaders(expression.getRowHeaders());
		MatrixStruct keepGenes = new MatrixStruct(expression.rows(),1);
		for(int r =0; r < expression.rows(); r++)
		{
			double gM = expression.sumRow(r)/expression.cols();
			if(Double.isFinite(gM))
			{
				geoMean.matrix.set(r,0,gM);
				geoMean.matrix.set(r,0,Math.pow(10,geoMean.matrix.get(r,0)));
				keepGenes.setRowHeader(r, expression.getRowHeaders()[r]);
			}
		}
		
		geoMean.keepRows(keepGenes);
		if(writeGeoFN == null)
			writeGeoFN = writeFolder+ "geoMean.txt";
		p("geofilename=" + writeGeoFN);
		this.geoFN=writeGeoFN;
		geoMean.write(writeGeoFN);
		
		return geoMean;
	}
	
	RlogLargeMatrix checkArgs(String[] args) 
	{
		if(args.length < 1)
		{
			p("Script requires the following arguments:\n"
					+ "1. filename=<expressionFN.txt> - Expression file with genes on rows samples on columns\n"
					+ "2. writeFolder=<writeFolderFN.txt> - Folder where the files will be written (default=parentFolder(input.txt))\n"
					+ "3. geoFN=<geoFn.txt> - Optional file with geometric mean per gene to use (default=null)\n");
			writeVars();
			System.exit(1);
		}
		
		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase())
			{
				case "json":
					;
				break;
				case "filename":
					expressionFN =value;
					break;
				case "writefolder":
					writeFolder = value;
					break;
				case "geofn":
					geoFN = value;
					break;
				case "writeall":
					writeAll = Boolean.parseBoolean(value);
					break;
				case "totalreadcount":
					logAdd = Double.parseDouble(value);
				break;
				case "roundvalues":
					roundValues = Boolean.parseBoolean(value);
				break;
				default:
					p("Incorrect argument supplied:\n"+ args[a] +"\nexiting");
					System.exit(1);
			}
		}
		writeVars();
		return this;
	}
}
