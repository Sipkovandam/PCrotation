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
	//A class that holds all the variables for PCcorrection and CreateGeneEigenvectorFile
	public  String jsonFN = "E:/Groningen/Data/Juha/Genes31995/Healthy/PCA/config.json";
	public  String expFile = "/Volumes/Promise_RAID/GeneNetwork/Sipko/TestExpression2/TESTexpression.txt";//expression file
	public  String chromLocationsFile = null;//File that contains the chromosome locations
	public  String writeFolder = "E:/Groningen/Data/Juha/Genes31995/31.07.pc1.illumina.genes.expressed.DEseqnorm/PCcorrection/";//Folder where to write
	public  String xmlFN= "config.xml";//don't touch this
	public  String GCgenes = null;//Not tested, don't use. If a file is supplied with the GC content per gene it will correct for this 
	public  String removeGene = null;//if there is a particular gene you want to remove from the matrix for any reason
	public  String genesToInclude = null;//List of genes. This list is selected prior to any other steps. Other genes are discarded
	
	public  boolean correlation = false; //if false uses covariance
	public  boolean setLowestToAverage = false;// sets all the lowest values in a sample to the average, effectively this means any gene that has an expression of does not contribute to the covariance or correlation
	public  boolean centerSamples = false;// adjusts the expression of genes so that the average expression of a sample become 0 (centering over the samples)
	public  boolean centerGenes = true;// adjusts the expression of genes so that the average expression of a sample become 0 (centering over the samples)
	public  boolean writeAll = true; //write all intermediate files, is slower but helps finding understanding what happens in each step
	public  boolean correctInputForSTdevs = false;//corrects the input for the standard deviation, can be done over genes or samples (look at the function itself it has a true/false argument)
	public  boolean correctInputForSTdevsAfterCenter = false;//same as previous, but this time after centering the data
	public  boolean log2 = true;
	public  boolean skipQuantileNorm = true;	
	public  boolean STdevCutoff = false;//calculates the standard deviation of all genes and throws out the <highestExpressed> percent genes with the lowest standard deviation (instead of using average expression for cutoff)
	public  boolean zScores = false;
	public  boolean directPCA = true;//PCA over genes rather than samples.
	public  boolean rLog = true;//If true uses DEseq normalization
	
	public  int minSamplesExpressed = -1;// if -1, this will include all samples. Otherwise will exclude any gene that is expressed in less then this number of samples (and those that have no variance)
										 // any gene that has an expression lower then "minExpression" is considered as "not expressed"
	public  int minExpression = -1; //if -1 does nothing. Otherwise sets the cutoff for minSamplesExpressed; see "minSamplesExpressed"
//	public  int ignoreLowestValues = -1;//Sets lowest values to the average so it does not contribute toward positive correlation.
	
	public  double addLogVal = 0.5; //Value to add before taking the logarithm
	public  double correctTotalReadCount = 0;//log((gene+0.5)/total*value)

	public  double randomValue = 0;//if 0 does nothing. Otherwise adds a random value that is below this value, to values below this value.
	public  double duplicateCutoff = 1; //if 1 does nothing. Samples with a correlation above this value are removed (only samples next to each other are removed, as this is usually the case for replicates, which we aim to remove like this, but not others)
	public  double highestExpressed = 1;//1 = all genes, 0.5 = 50% highest expressed genes only (removes 50% lowest expressed genes after quantile normalization (then re-normalizes)).
	public  double spearman = -1;//if 0, does spearman, if >0, it sets all values below this value to 0.
	
	public  String correlationScript = "/Volumes/Promise_RAID/GeneNetwork/Sipko/CorrelationLargeGenes.jar"; //script used to calculate correlation.
	
	//Prediction script
	public  String geneNameFile = "/Volumes/Promise_RAID/GeneNetwork/ENSGToGeneNameHGNCBiotypeChromosomeStartStopStrandAndDescriptionV75.txt.filtered.txt";
	public  String itemSetFile = "/Volumes/Promise_RAID/GeneNetwork/Sipko/Scripts/HPO/HPO.gmt";
	public  String geneTermOutFile = "HPO_geneterm.txt";
	public  String termDesCoutFile = "HPO_terms.txt";
	public  String itemType = "Gene";
	public  String limitToItemsInGenesetFile = "false";//Data/OldSchoolExpression/GeneAnnotation/Genes.txt
	public  String maxItems = "10000";//Maximum number of items allowed in Tessa script
	public  String minItems = "10";//Not really sure how this is different from absMinItems
	public  String absMinItems = "5";//absolute minimum number of items in geneterm
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
	public  String avgStdevFolder = null;//this is the folder with the files for the z-score calculations containing the averages and standard deviations to be used for this.
										 //allows z-scores to be calculated using the BBMRI averages and stdevs
	public  double zScoresCutoff = Double.parseDouble("0");
	public  boolean correctResultsForSTdevs = true;
	
	public int optimalPCremoval = -1;
	public String tempName;
	
	
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
	
}
