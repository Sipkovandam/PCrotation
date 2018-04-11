package GeneNetwork;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.math3.stat.inference.WilcoxonSignedRankTest;
import org.biojava.ontology.Ontology;
import org.biojava.ontology.Term;
import org.biojava.ontology.Triple;
import org.biojava.ontology.io.OboParser;

import com.sun.msv.datatype.xsd.datetime.ParseException;

import MatrixScripts.MyMatrix;
import Tools.FileUtils;
import Tools.Script;
import umcg.genetica.math.stats.WilcoxonMannWhitney;

public class ChildHpoBasedZscorePredictor extends Script<ChildHpoBasedZscorePredictor>
{
	//create a z-score matrix based on weighted child z-scores for each HPO term
	
	//starts in the top of the tree and works its way down, and then back up recusively calculating z-scores.
	//only uses genes that are present in HPO terms as background
	//when calculating weighted z-scores based on a child, only uses children with significant AUCs (corrected for bonferroni)
	//
	
	transient Ontology hpoOntology;
	transient String[] genes = null;
	transient HashMap<String, HashMap<String, Integer>> hpoToGenes= null; 
	transient HashMap<String, String> ensembl_To_GeneSymbol= null; 
	transient HashMap<String, String> geneSymbol_To_Ensembl= null; 

	transient MyMatrix geneToZscore = null;
	transient BufferedWriter aucWriter = null;
	transient BufferedWriter newZscoreWriter = null;
	transient HashMap<String, Double[]> hpo_To_CalculatedChildZscores = new HashMap<String, Double[]>();
	
	String hpoTreeFn="";
	String hpoToGeneFn="";
	String geneHpoZscoreMatrixFn="";
	String ensembl_To_GeneSymbolFn="";

	String startHpo="HP:0000118";
	
	String newZscoreWriteFn = null;
	
	boolean onlyUseChildGenes = true;
	
	@Override
	public void run()
	{
		try
		{
			init();

			Term startTerm = hpoOntology.getTerm(startHpo);			
			Term is_a = hpoOntology.getTerm("is_a");
			Double[] weightedZscore = childBasedAnalysis(startTerm, is_a,1);
			
			//save new zscore file file
			aucWriter.close();
			newZscoreWriter.close();
		}catch(Exception e){e.printStackTrace();}
		
		
	}

	private void init() throws FileNotFoundException, IOException, org.biojava.bio.seq.io.ParseException, ParseException
	{
		log("Reading ensembl to gene symbol hash");
		ensembl_To_GeneSymbol= FileUtils.readStringStringHash(ensembl_To_GeneSymbolFn);
		geneSymbol_To_Ensembl= FileUtils.readStringStringHash(ensembl_To_GeneSymbolFn,1,0);
		
		log("Reading hpoToGene matrix");
		hpoToGenes=FileUtils.readStringMultiStringHash(hpoToGeneFn, 0, 3, false);
		log(hpoToGenes.size() + " hpo terms loaded");

		//read the z-score matrix;
		log("Reading geneToZscore matrix:\t" + geneHpoZscoreMatrixFn);
		geneToZscore= new MyMatrix(geneHpoZscoreMatrixFn);
		
		log("Create fileWriters");
		newZscoreWriter=FileUtils.createWriter(newZscoreWriteFn);
		newZscoreWriter.write(FileUtils.StringArrayToWriteString(geneToZscore.colNames)+"\n");
		aucWriter=FileUtils.createWriter(FileUtils.removeExtention(newZscoreWriteFn)+"_AUCs.txt");
		aucWriter.write("Hpo term\told AUC\told p-value\tnew AUC\tnew p-value\n");
						
		//get all genes in HPO terms
		log("Reading genes to include");
		HashSet<String> genesHashSet = getAllGenes(hpoToGeneFn);
		genes = genesHashSet.toArray(new String[genesHashSet.size()]);
		
		log("Loading hpo tree:\t" + hpoTreeFn);
		hpoOntology = loadHpoOntology(new File(hpoTreeFn));
		
	}

	private HashSet<String> getAllGenes(String phenotypeToGeneFn) throws FileNotFoundException, IOException
	{
		BufferedReader phenotypeToGeneReader = FileUtils.createReader(phenotypeToGeneFn);
		String line = phenotypeToGeneReader.readLine();//skip header
		HashSet<String> genes = new HashSet<String>();
		while((line = phenotypeToGeneReader.readLine())!=null)
		{
			genes.add(line.split("\t")[3]);
		}
		
		return genes;
	}

	private Double[] childBasedAnalysis(Term hpoTerm, Term is_a, int hpoTreeLevel) throws IOException
	{
		//find all child terms
		Set<Triple> childHpoTerms =hpoOntology.getTriples(null, hpoTerm, is_a);
		//this is where the z-scores of the child terms will be stored
		ArrayList<Double[]> subHpoTermZscoresArrayList= new ArrayList<Double[]>();
		
		HashSet<String> genesInHpo = new HashSet<String> ();
		for(Triple child: childHpoTerms)
		{
			log(child.getName() + "\thpoTreeLevel=\t" + hpoTreeLevel + "\thpoTerm:\t" + hpoTerm.getName());
			hpoTreeLevel++;
			
			//for each child determine the z-scores based on their own children
			Double[] subHpoTermZscores=hpo_To_CalculatedChildZscores.get(hpoTerm.getName());
			if(!hpo_To_CalculatedChildZscores.containsKey(hpoTerm.getName()))
				subHpoTermZscores=childBasedAnalysis(child.getSubject(), is_a, hpoTreeLevel);
			
			if(subHpoTermZscores!=null)
			{
				subHpoTermZscoresArrayList.add(subHpoTermZscores);
				log("child1=" + child.getSubject().getName());
				log("genes=" + hpoToGenes.get(child.getSubject().getName()));
				if(hpoToGenes.get(child.getSubject().getName())!=null)
				{
					log("child=" + child.getSubject().getName() + "\t" + hpoToGenes.get(child.getSubject().getName()));
					genesInHpo.addAll(hpoToGenes.get(child.getSubject().getName()).keySet());
				}
			}
		}
		
		log("Calculating new z-scores for hpoTerm:\t " + hpoTerm.getName());
		
		log("hpoTerm.getName()=" + hpoTerm.getName() + "\t" + hpoToGenes.get(hpoTerm.getName()));
		if(genesInHpo.size() == 0 && hpoToGenes.get(hpoTerm.getName())!=null)
			genesInHpo.addAll(hpoToGenes.get(hpoTerm.getName()).keySet());
		
//		if(hpoToGenes.get(hpoTerm.getName())!=null)
//			genesInHpo.addAll(hpoToGenes.get(hpoTerm.getName()).keySet());

		Double[] childBasedZscores= null;
		double[] newAuc = new double[]{0.5,1};	
		//if it has children (with significant AUCs) calculate new weighted z-scores based on the children
		if(subHpoTermZscoresArrayList.size()>1)	//could make this 1...	
		{
			childBasedZscores=calculateWeightedZscore(subHpoTermZscoresArrayList);
			//calculate the AUC for the new weighted z-scores
			HashMap<String, Integer> genesInterm = hpoToGenes.get(hpoTerm.getName());
			if(onlyUseChildGenes)
			{
				if(genesInHpo!=null && genesInterm!=null)
					log("Genes in children= " + genesInHpo.size() + "\tgenes in parent= " + genesInterm.size());
				genesInterm=geneTermRemoveNonChildGenes(genesInHpo);
				if(genesInHpo!=null && genesInterm!=null)
					log("genesInterm= " + genesInterm.size());
			}
				
			newAuc = getAucForHpoTerm(childBasedZscores, genesInterm);
		}
		
		//get the old z-scores
		Double[] zScores= geneToZscore.getColValues(hpoTerm.getName(), true);
		//calculate HPO based on others
		double[] oldAuc = getAucForHpoTerm(zScores, hpoToGenes.get(hpoTerm.getName()));

		//write the auc scores for the old and the new calculations
		int size=0;
		if( hpoToGenes.get(hpoTerm.getName())!=null)
			size= hpoToGenes.get(hpoTerm.getName()).size();
		String aucLine = hpoTerm.getName()+"\t" + oldAuc[0] + "\t" +oldAuc[1] + "\t" + newAuc[0]+ "\t" + newAuc[1]+"\t" +size;
		aucWriter.write(aucLine + "\n");
		log("AucLine=\t" + aucLine);
		
		//determine the new z-scores to return for this child term
		double bonferoniCutoff=0.05/(double)geneToZscore.cols();
		Double[] newZscores= zScores;
		if(oldAuc[0] < newAuc[0] && newAuc[1]<bonferoniCutoff)
			newZscores=childBasedZscores;
		
		if(newZscores!=null)
			newZscoreWriter.write(hpoTerm.getName()+FileUtils.doubleArrayToWriteString(newZscores)+"\n");
		
		//if neither the new or the old AUC is signicant do not use this term to calculate the weighted Z-scores for the parents.
		if(oldAuc[1] > bonferoniCutoff && newAuc[1] > bonferoniCutoff)
			return null;
		
		hpo_To_CalculatedChildZscores.put(hpoTerm.getName(), newZscores);
		return newZscores;
	}

	private HashMap<String, Integer> geneTermRemoveNonChildGenes(HashSet<String> genesInChildren)
	{
		HashMap<String, Integer> genesInChildrenGeneSymbol =new HashMap<String, Integer>();
		for(String gene: genesInChildren)
				genesInChildrenGeneSymbol.put((gene), 0);
		
		return genesInChildrenGeneSymbol;
	}

	private double[] getAucForHpoTerm(	Double[] zScores,
							HashMap<String, Integer> hpoGenes)
	{
		
		
		if(hpoGenes==null || hpoGenes.size()<3)
			return new double[]{0.5,1};
		
		int nHpo =0;
		int nOthers = 0;
		//get the array sizes needed (i think it is faster than using arraylist<Double>)
		for(String gene : geneToZscore.rowNames)
		{
			if(ensembl_To_GeneSymbol.get(gene)!=null && hpoGenes.containsKey(ensembl_To_GeneSymbol.get(gene)))
			{
				nHpo++;
			}
			else
			{
				nOthers++;
			}
		}
		
		double[] zScoresHpoGenes= new double[nHpo];
		double[] zScoresOtherGenes= new double[nOthers];

		int x = 0;
		int xHpo =0;
		int xOthers = 0;
		for(String gene : geneToZscore.rowNames)
		{
			if(hpoGenes.containsKey(ensembl_To_GeneSymbol.get(gene)))
			{
				zScoresHpoGenes[xHpo]=zScores[x];
				xHpo++;
			}
			else
			{
				zScoresOtherGenes[xOthers]=zScores[x];
				xOthers++;
			}
			x++;
		}
		
		double[] auc = new double[]{0.5,1};
		WilcoxonMannWhitney wilcoxonMannWhitney = new WilcoxonMannWhitney();
		if(zScoresHpoGenes.length>1)
		{
			auc[1] = wilcoxonMannWhitney.returnWilcoxonMannWhitneyPValue(zScoresHpoGenes, zScoresOtherGenes);
			auc[0] = wilcoxonMannWhitney.getAUC();
		}
			
		return auc;
	}

	private Double[] calculateWeightedZscore(ArrayList<Double[]> subHpoTermZscores)
	{
		Double[] totalPerGene = new Double[subHpoTermZscores.get(0).length];
		for(Double[] zScores : subHpoTermZscores)
			for(int z=0;z < zScores.length; z++)
			{
				if(totalPerGene[z]==null)
					totalPerGene[z]=(double) 0;
				totalPerGene[z]+=zScores[z];
			}
		
		double divider =Math.sqrt(subHpoTermZscores.size());
		for(int z=0;z < totalPerGene.length; z++)
			totalPerGene[z]+=totalPerGene[z]/divider;
		
		return totalPerGene;
	}

	public Ontology loadHpoOntology(final File hpoOboFile) throws FileNotFoundException, IOException, ParseException, org.biojava.bio.seq.io.ParseException {
		OboParser parser = new OboParser();
		BufferedReader oboFileReader = new BufferedReader(new InputStreamReader(new FileInputStream(hpoOboFile)));
		Ontology hpoOntology = parser.parseOBO(oboFileReader, "HPO", "HPO");
		return hpoOntology;
	}

}

