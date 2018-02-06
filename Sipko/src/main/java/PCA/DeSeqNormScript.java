package PCA;

import java.io.File;
import java.io.IOException;

import MatrixScripts.MyMatrix;
import Tools.FileUtils;
import Tools.Script;

public class DeSeqNormScript extends Script<DeSeqNormScript>
{
	//This is doing the DESeq normalization, which does not use replicate information.
	
	String expressionFN = "E:/Groningen/Scripts/Tests/Rlog.java/Samples.txt";
	
	String writeFolder = null;//if null becomes new File(expressionFN).getParent()+"/";
	String geoFN = null;//if null calculates geometric means based on this dataset
	boolean writeAll = true;	//write all intemediary files too
	double logAdd = 1;//previously used to multiply the results by this number, but seems pointless since it does not have any effect, so does not do anything anymore
	boolean log = false;
	boolean genes = true;
	boolean roundValues = true;//rounds expression values to whole counts
	
	@Override
	public void run()
	{
		try
		{
			init();
		
			DeSeqNorm.main(new String[3]);
		} catch (IOException e){e.printStackTrace();}
	}

	private void init()
	{
		DeSeqNorm.expressionFN=this.expressionFN;
		DeSeqNorm.writeFolder=this.writeFolder;
		DeSeqNorm.geoFN=this.geoFN;
		DeSeqNorm.writeAll=this.writeAll;
		DeSeqNorm.logAdd=this.logAdd;
		DeSeqNorm.log=this.log;
		DeSeqNorm.genes=this.genes;
		DeSeqNorm.roundValues=this.roundValues;	
		if(writeFolder == null)
			writeFolder = new File(expressionFN).getParent()+"/DESeqNorm/";
	}
}
