package GeneNetwork;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;

import org.apache.commons.io.FileUtils;

import Tools.Script;

public class ExomizerAnalysis extends Script<ExomizerAnalysis>
{
	//creates config files for exomizer analysis for all vcf files in "inputFolder"
	
	String dna_To_HpoFn = "/groups/umcg-gdio/tmp03/umcg-svandam/Data/VCF/Mixed/Dna_To_HPO.txt";
	String inputFolder = "/groups/umcg-gdio/tmp03/umcg-svandam/Data/VCF/Mixed/";
	String outputFolder = "/groups/umcg-gdio/tmp03/umcg-svandam/Data/VCF/Mixed/ExomizerSettings/";
	String templateFn = "/groups/umcg-gdio/tmp03/umcg-svandam/Data/VCF/Mixed/ExomizerSettings/analysis-exome-template.yml";
	
	public void run()
	{
		try
		{
			String templateFile=FileUtils.readFileToString(new File(templateFn));
			HashMap<String, String> dna_To_Hpo = Tools.FileUtils.readStringStringHash(dna_To_HpoFn);
			log("dna_To_Hpo =\t" + dna_To_Hpo);
			File[] files = new File(inputFolder).listFiles(); 
			for(File file : files)
			{
				if(!file.getName().contains(".vcf"))
					continue;
				
				String fn = file.getName();
				String dnaNumber= fn.split("_noInfo")[0];
				log("DNA number =\t" + dnaNumber);
				String hpoTerms = dna_To_Hpo.get(dnaNumber);
				if(hpoTerms==null)
				{
					log("No hpoTerms available for this sample:" + dnaNumber);
					continue;
				}
				
				String outputFile = templateFile.replace("vcf: fileName ", "vcf: "+ file.getAbsolutePath());
				String outputDir=Tools.FileUtils.makeFolderNameEndWithSlash(outputFolder)+dnaNumber+"/";
				Tools.FileUtils.makeDir(outputDir);
				outputFile=outputFile.replace("hpoIds: hpoIds", "hpoIds: "+ hpoTerms);
				outputFile=outputFile.replace("results/Pfeiffer-hiphive-exome-PASS_ONLY", outputDir);
				
				String writeFn=outputFolder+dnaNumber+".yml";
				FileUtils.write(new File(writeFn), outputFile);
			}
			
		} catch (Exception e)
		{
			e.printStackTrace();
		}
		log("Finished, files written at:\t" + outputFolder);
	}
}
