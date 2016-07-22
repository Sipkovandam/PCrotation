package PCA;

import java.io.IOException;

import org.apache.commons.math3.stat.inference.TTest;

import Analyses.WilcoxonMannWhitney;
import JuhaPCA.PCA;

public class IdentifySimilarSamples 
{
	public static void main (String[] args) throws IOException
	{
		//The samples you want to rank (corrected samples)
//		String samplesFN = "E:/Groningen/Data/PublicSamples/Test10/10000Samples_Duplicate1.0_RndVal2_0.20Highest_log2/18DownSyndrome26Normal2Cancer_countsAdj/"
//				+ "PC_1-1000_DevidedBySTdevs.txt";
//		
//		//fileName that contains the genes in the rownames. These genes are the genes you expect to have an abberant expression (e.g. chr21 in down syndrome)
//		String genesFN   = "E:/Groningen/Data/PublicSamples/Test10/10000Samples_Duplicate1.0_RndVal2_0.20Highest_log2/18DownSyndrome26Normal2Cancer_countsAdj/"
//				+ "PC_1-1000_DevidedBySTdevsdiseaseGenes_19samples.txt";
		//E:\Groningen\Data\Iris\CountsSD200\directPCA_Rlog_0.2_covar\est_counts_nocancernocellline\
		//String samplesFN = "E:/Groningen/Data/Iris/CountsSD200/directPCA_Rlog_0.2_covar/est_counts_nocancernocellline/PC_1-300_DevidedBySTdevs.txt";
		String samplesFN = "E:/Groningen/Data/PublicSamples/Test13/directPCA_Voom_0.2/18DownSyndrome26Normal2Cancer_counts/PC_1-300_DevidedBySTdevs.txt";
		
		//String genesFN = "E:/Groningen/Data/Monogenetic_disease_samples_RadBoud/directPCA_Voom_0.2/CountsGENES_Radboud/PC_1-300_DevidedBySTdevsdiseaseGenes_7samples_3.0STdevs.txt";
		String genesFN = "E:/Groningen/Data/PublicSamples/Test13/directPCA_Voom_0.2/18DownSyndrome26Normal2Cancer_counts/19-04-2016/PC_1-300_OutliersCol0-2_3Stdevs.txt";
		String writeFN = "E:/Groningen/Data/PublicSamples/Test13/directPCA_Voom_0.2/18DownSyndrome26Normal2Cancer_counts/19-04-2016/PC_1-300_OutliersCol2,13,18_3Stdevs_Similar_NoChr21.txt";
		String removeGenesFN = "E:/Groningen/Data/GenePositionInfo_Chr21Only.txt";
		
		if(args.length<2)
			checkArgs(args);
		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase()){
				case "filename":
					samplesFN = value;
					break;
				case "genes":
					genesFN = value;
					break;
				case "writefile":
					writeFN = value;
					break;
				default:
					checkArgs(args);
					System.out.println("Incorrect argument supplied; exiting");
					System.exit(1);
			}
		}
		if(writeFN == null)
			writeFN = samplesFN.replace(".txt", "_")+ samplesFN.split("/")[samplesFN.split("/").length-2] + "_diseaseSimilar.txt";
		System.out.println("WriteFN = " + writeFN);
		String folderName = samplesFN.substring(0, samplesFN.lastIndexOf("/")+1);
		MatrixStruct samples = new MatrixStruct(samplesFN);
		if(samples.getColHeaders()[0].contains("ENSG000") || samples.getColHeaders()[0].contains("ENST000"))
		{
			System.out.println("Genes should be on rows, Transposing");
			samples.transpose();
		}
			
		MatrixStruct genes = new MatrixStruct(genesFN);
		MatrixStruct removeGenes = new MatrixStruct(removeGenesFN);
		genes.removeRows(removeGenes);
		genes.write(genesFN.replace(".txt", "GENESKEPT.txt"));
		
		MatrixStruct output = new MatrixStruct(samples.cols(), 8);
		output.setRowHeaders(samples.getColHeaders());
		output.setColHeaders(new String[]{"pValue","AverageDisease","AverageOthers","cutoff","STdev", "dups+dels","wilcoxon p", "wilcoxon AUC"});
		
		String getSamples = null;
		MatrixStruct smoothed = new MatrixStruct(samplesFN.substring(0, samplesFN.lastIndexOf("_"))+"_.Smoothed101Genes.txt");
		
		for(int s = 0; s < samples.cols(); s++)
		{
			int d = 0;
			int o = 0;
			double[] disease = new double[genes.rows()];
			double[] others = new double[samples.rows()-genes.rows()];
			double[] all = new double[samples.rows()];
			for(int r = 0; r < samples.rows(); r++)
			{
				if(genes.rowHash.containsKey(samples.getRowHeaders()[r]))
				{
					int row = genes.rowHash.get(samples.getRowHeaders()[r]);
					double diseaseValue = genes.matrix.get(row, 0);
					if((samples.matrix.get(r, s)< 0 && diseaseValue > 0) || (samples.matrix.get(r, s)> 0 && diseaseValue < 0))
						disease[d] = samples.matrix.get(r, s);//punishment for being opposite ;)
					else
						disease[d] = Math.abs(samples.matrix.get(r, s));
					d++;
				}
				else
				{
					others[o] = Math.abs(samples.matrix.get(r, s));
					o++;
				}
				all[r] =Math.abs(samples.matrix.get(r, s));
			}
			double diseaseAverage = org.apache.commons.math3.stat.StatUtils.mean(disease);
			double othersAverage = org.apache.commons.math3.stat.StatUtils.mean(others);
			double stdev =  Math.pow(org.apache.commons.math3.stat.StatUtils.variance(all), 0.5);
			double tTest = new TTest().tTest(disease,others)/2;
			if(othersAverage > diseaseAverage)//This tends to happen in samples with a lot of deletions and can lead to very low p-values
				tTest=1;
			output.matrix.set(s, 0, tTest);
			output.matrix.set(s, 1, diseaseAverage);
			output.matrix.set(s, 2, othersAverage);
			output.matrix.set(s, 3, 0.20);
			output.matrix.set(s, 4, stdev);
			WilcoxonMannWhitney wmw = new WilcoxonMannWhitney(); 
			double pValue = wmw.returnWilcoxonMannWhitneyPValue(disease, others); 
//			System.out.println("dis");
//			for(int dis = 0 ;dis<disease.length;dis++)
//				System.out.println(disease[dis]);
//			System.out.println("Others");
//			for(int dis = 0 ;dis<others.length;dis++)
//				System.out.println(others[dis]);
			
			double auc = wmw.getAUC();
			//System.out.println("mann pval =" + pValue + " auc=" + auc);
			output.matrix.set(s, 6, pValue);
			output.matrix.set(s, 7, auc);
			
			countDuplications(output, s, smoothed);
		}
		output.sortCol(0,-1);
		output.write(writeFN);
		PCA.log("Done, file written to: " + writeFN);
		samples = null;
		
		int n = 0;
		for(int r = 0; r < output.rows(); r++)
			if(output.matrix.get(r, 0) < 0.25)
			{
				if(getSamples == null)
					getSamples = output.getRowHeaders()[r];
				else
					getSamples += ","+output.getRowHeaders()[r];
			}

			
		System.out.println(samplesFN);
//		GetCols.main(new String[]{"fileName="+samplesFN,"getgenes="+getSamples,"writename="+samplesFN.replace(".txt", "_diseaseSamples"+getSamples.split(",").length+".txt")});
//		GetCols.main(new String[]{"fileName="+folderName+"PC_1-0_DevidedBySTdevs.txt","getgenes="+getSamples,"writename="+samplesFN.replace(".txt", "PC_1-0_DevidedBySTdevs_diseaseSamples"+getSamples.split(",").length+".txt")});
//		GetCols.main(new String[]{"fileName="+folderName+"PC_1-300_.Smoothed101Genes.txt","getgenes="+getSamples,"writename="+samplesFN.replace(".txt", "PC_1-300_.Smoothed101Genes_diseaseSamples"+getSamples.split(",").length+".txt")});
	}

	private static void countDuplications(MatrixStruct output, int s, MatrixStruct smoothed) 
	{
		int n = 0;
		double stdev = output.matrix.get(s, 3);
		int gap = 0;
		for(int g = 0; g < smoothed.rows(); g++)
		{
			double value = smoothed.matrix.get(g, s);	
			//basically if the average of 25 genes is > 1 stdev away the we consider this region as duplicated (or deleted)
			if((value > stdev || value < -stdev))//50 is minimum gap between genes before it is considered another cytoband 
													 //(so 50 genes have to be below stdev)(also need at least 50 genes above stdev before considered a duplication (or deletion))
			{
				if(gap >= 50)
				{
					n++;	
				}
				gap = 0;
			}
			else
				gap++;
		}
		output.matrix.set(s,5,n);
	}

	private static void checkArgs(String[] args) {
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
			return;
		System.out.println("This script identifies outlier genes across a number of samples\n"
				+ "fileName=<fileName> - FileName containing the corrected samples\n"
				+ "genes=<savefileName> - Name of the file containing the outlier/disease genes\n"
				+ "writefile=<number> - Filename to save the files to (optional)\n");
		System.exit(1);
	}
}
