package PCA;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.w3c.dom.Document;
import org.w3c.dom.Element;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

public class Var 
{
	public  boolean correlation = true; //if false uses covariance
	public  boolean setLowestToAverage = false;// sets all the lowest values in a sample to the average, effectively this means any gene that has an expression of does not contribute to the covariance or correlation
	public  boolean adjustSampleAverages = true;// addjusts the expression of genes so that the average expression of a sample become 0 (centering over the samples)
	public  boolean writeAll = true; //write all intermediate files, is slower but helps finding understanding what happens in each step
	public  boolean correctInputForSTdevs = false;//corrects the input for the standard deviation, can be done over genes or samples (look at the function itself it has a true/false argument)
	public  boolean correctInputForSTdevsAfterCenter = false;//same as previous, but this time after centering the data
	public  boolean log2 = true;
	public  boolean skipQuantileNorm = true;	
	public  boolean STdevCutoff = false;
	public  boolean zScores = false;
	public  boolean directPCA = true;
	
	public  int minSamplesExpressed = -1;// if left -1 and minExpression>0, this will become all samples
	public  int minExpression = -1; //if left -1 does nothing
	public  int ignoreLowestValues = -1;//I still need to fix this
	
	public  double addLogVal = 1;
	public  double correctTotalReadCount = 0;//log((gene+0.5)/total*value)
	public  double rLog = 1000000;
	public  double topVariance = 1;
	public  double randomValue = 0;
	public  double duplicateCutoff = 1;
	public  double highestExpressed = 1;//1 = all genes, 0.5 = 50% highest expressed genes only (removes 50% lowest expressed genes after quantile normalization (then re-normalizes)).
	public  double spearman = -1;//if 0, does spearman, if >0, it sets all values below this value to 0.
	
	public  String jsonFN = "config.json";
	public  String removeGene = null;
	public  String expFile = null;
	public  String chromLocationsFile = null;
	public  String writeFolder = null;
	public  String xmlFN= "config.xml";
	public  String correctGCSamples = null;
	public  String GCgenes = "";

	public  String correlationScript = "/Volumes/Promise_RAID/GeneNetwork/Sipko/CorrelationLargeGenes.jar";
	
	//Prediction script
	public  String geneNameFile = "/Volumes/Promise_RAID/GeneNetwork/ENSGToGeneNameHGNCBiotypeChromosomeStartStopStrandAndDescriptionV75.txt.filtered.txt";
	public  String itemSetFile = "/Volumes/Promise_RAID/GeneNetwork/Sipko/Scripts/HPO/HPO.gmt";
	public  String geneTermOutFile = "HPO_geneterm.txt";
	public  String termDesCoutFile = "HPO_terms.txt";
	public  String itemType = "Gene";
	public  String limitToItemsInGenesetFile = "false";//Data/OldSchoolExpression/GeneAnnotation/Genes.txt
	public  String maxItems = "10000";
	public  String minItems = "10";
	public  String absMinItems = "5";
	public  String useTtest = "true";//if false Wilcoxon will be used
	public  String runRealAnalysis = "true";
	public  String nrPermutations = "0";
	public  String writeBinary = "false";
	public  String logToFile = "false";
	public  String label = "HPO-unscaled";
	public  String javaExe = "/usr/bin/java -jar";
	public  String RNAseqJar = "/Volumes/Promise_RAID/juha/PublicRNASeq/Prediction/RNASeq.jar";
	public  String populateGenesetDBjs = "/Volumes/Promise_RAID/GeneNetwork/Sipko/05-2016/Tessa/populateGenesetDBTXT.js";
	public  String geneDB = "genedb";
	public  String hpoDB = "hpodb";
	public  String termtype = "HPO";
	public  String tessaPredictionFolder = "/Volumes/Promise_RAID/GeneNetwork/Sipko/05-2016/Tessa/Server_testing_DB/";
	public  String tessaPredictionScript = "Ranking_multipleVersionsCombined.js";
	public  String tessaMaxGenesPerTerm = "max15";
	public  String tessaType = "regular";
	
	//if these variables are set the PC correction will be run after the eigenvectors are created (or if eigenvectors already exist only this part is run)
	public  String sampleFile = "";
	public  String chr21FN = "";
	public 	String PCs = "";
	public  String writeFolderCorrected = "";
	public  double zScoresCutoff = Double.parseDouble("0");
	public  boolean correctResultsForSTdevs = true;

	public int optimalPCremoval = -1;
	
	public void writeVars()
	{
		writeVars(this.jsonFN);
	}
	public void writeVars(String jsonFN)
	{
		this.jsonFN=jsonFN;
		Gson gson = new GsonBuilder().setPrettyPrinting().create();
		//System.out.println(gson.toJson(this));
		//System.out.println(gson.toString());
		System.out.println(jsonFN);
		System.out.println("writefolder = " + writeFolder);
		System.out.println(getWritePath(jsonFN));
		
		try
		{
			BufferedWriter writer = new BufferedWriter(new FileWriter(new File(getWritePath(jsonFN))));
			
			writer.write(gson.toJson(this));
			writer.close();
		}catch(Exception e){e.printStackTrace();}
	}
	public void readVars()
	{
		readVars(this.jsonFN);
	}
	public Var readVars(String jsonFN)
	{
		String jsonString = "";
		
		try
		{
			File file = new File(jsonFN);
			if(!file.exists())
			{
				System.out.println("Json file does not exist:" + jsonFN);
				System.exit(0);
			}
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String line = null;
			while((line=reader.readLine())!=null)
			{
				jsonString+=line;
			}
			reader.close();
		}catch(Exception e){e.printStackTrace();}
		
		Gson gson = new GsonBuilder().setPrettyPrinting().create();
		Var var = gson.fromJson(jsonString, Var.class);
		System.out.println(gson.toJson(var));
		return var;
	}

	//Prediction script
	public boolean filePathsExist()
	{
		File[] files = new File[7]; 
				
		files[0] =new File(this.expFile);
		files[1] = new File(this.geneNameFile);
		files[2] = new File(this.itemSetFile);
		files[3] = new File(this.RNAseqJar);
		files[4] = new File(this.populateGenesetDBjs);
		files[5] = new File(this.tessaPredictionFolder);
		files[6] = new File(this.correlationScript);
		
		for(File file : files)
		{
			if(file.exists())
				continue;
			System.out.println("THIS FILE/FOLDER DOES NOT EXIST!\n" + file.getAbsolutePath());
			return false;
		}
		return true;
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
}
