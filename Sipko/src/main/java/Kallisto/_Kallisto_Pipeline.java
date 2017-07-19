package Kallisto;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import MatrixScripts.LaneMerger;
import STAR.STAR_Variables;
import Slurm.Slurm;
import Tools.FileUtils;
import Tools.JSONutil;
import Tools.Script;
import Tools.FileSearcher;

public class _Kallisto_Pipeline extends Script<_Kallisto_Pipeline> {
	// This class does the following:
	// 1. Finds all fastQ files in a folder and subfolders
	// 2. Determines which fastQ files are pairs
	// 3. Maps the files using Kallisto
	// 4. Creates a file describing the number of reads per sample mapped
	// 5. Sums reads of transcritps to genes
	// 6. Creates 1 large matrix for genes and one for transcripts
	// 7. Creates 1 large matrix for samples >70% mapping

	String outputFolder = "E:/Groningen/Test/JSON/";
	String input = "J:/DATA/";//either 1. a comma separated list of folders from which to include all files or 2. a filename containing all the filenames of the files to be included
	String fastQSearchStrings = ".fastq,.fq";
	String forbiddenSearchStrings = ".md5";

	// Kallisto variables
	String kallistoOutputRoot = outputFolder + "Kallisto/";
	String kallistoVersion = "Kallisto/0.42.2.1-goolf-1.7.20";// "Kallisto/0.42.4-goolf-1.7.20";
	String kallistoThreads = "4";
	String kallistoWalltime = "05:59:00";
	String kallistoMaxMemory = "8gb";
	String pairStrings = "_R1,_R2,_1.fq,_2.fq";

	String slurmUserName = "umcg-svandam";
	String finishedemailaddress = "sipkovandam@gmail.com";

	// file merging Kallisto Files, creating expression per gene and
	// creating cutoff files variables
	int kallistoOutputColumn = 2;// 2 for counts, 3 for tpm values
	double kallistoThreshold = 0.7;
	String transcriptsToGenesFN = "";
	int minpercentagefeaturesexpressed = 0;
	private Kallisto_Variables kallisto_Variables = new Kallisto_Variables();
	
	private String sampleSheetFn = null;

	@Override
	public void run()
	{
		try
		{
			// args = new String[]{"json=E:/Groningen/Test/JSON/var.json"};
			String fastQFiles = input;
			if(input.contains(",") || !new File(input).isFile())
			{
				if(input.contains(".txt"))
				{
					p("THis file does not exist:\n" + input +"\n Exiting");
					System.exit(2);
				}
				fastQFiles = findFastQFiles();
				
			}
	
			kallistoSlurm(fastQFiles);
	
			String geneExpressionFn=combineKallistoOutput();
			
			//merge lanes
			if(sampleSheetFn!=null)
			{
				LaneMerger laneMerger = new LaneMerger();
				laneMerger.setCountsFn(geneExpressionFn);
				laneMerger.setSummedFn(FileUtils.removeExtention(geneExpressionFn) + "_LanesMerged.txt.gz");
				laneMerger.setLaneIndexInFileName(4);
				laneMerger.setSampleSheetFn(sampleSheetFn);
				laneMerger.run();
			}
			
			
		}catch(Exception e){e.printStackTrace();}
	}

	private String combineKallistoOutput() throws Exception 
	{
		String mappingPercentagesFN = kallistoOutputRoot + "mappingPerSample.txt";
		String tsvFilesToShScriptFN = kallistoOutputRoot + "scriptNumberToFiles.txt";
		File kallistoFolder = new File(kallistoOutputRoot);
		if(!kallistoFolder.exists())
			kallistoFolder.mkdirs();
		
		CombineKallisto combineKallisto = new CombineKallisto();
		combineKallisto.setGenesFN(this.kallistoOutputRoot+"counts_GENES.txt.gz");
		combineKallisto.setKallistoOutputFolder(this.kallistoOutputRoot);
		combineKallisto.setKallistoColumn(this.kallistoOutputColumn);
		combineKallisto.setTranscriptsToGenesFN(this.transcriptsToGenesFN);
		combineKallisto.setThreshold(this.kallistoThreshold);
		combineKallisto.setMinPercentageFeaturesExpressed(this.minpercentagefeaturesexpressed);
		combineKallisto.setMappingPercentagesFN(mappingPercentagesFN);
		combineKallisto.setTsvThatshouldBeThereFN(tsvFilesToShScriptFN);

		combineKallisto.run();
		
		return combineKallisto.getGenesFN();
	}


	private void kallistoSlurm(String fastQFiles) throws Exception {
		Slurm<Kallisto_Variables> kallistoSettings = new Slurm<Kallisto_Variables>(kallistoVersion);
		
		kallistoSettings.setFastQFiles(fastQFiles);
		kallistoSettings.setMapper(kallistoVersion);
		kallistoSettings.setOutputRoot(kallistoOutputRoot);
		kallistoSettings.setThreads(Integer.parseInt(kallistoThreads));
		kallistoSettings.setWalltime(kallistoWalltime);
		kallistoSettings.setMaxMemory(kallistoMaxMemory);
		kallistoSettings.setPairStrings(pairStrings.split(","));
		kallistoSettings.setSlurmUserName(slurmUserName);
		kallistoSettings.setFinishedemailaddress(finishedemailaddress);
		kallistoSettings.setMapperVars(this.kallisto_Variables);
		
		
		kallistoSettings.run();
	}

	private String findFastQFiles() throws Exception {
		String fastQFiles = outputFolder + "fastQfiles.txt";
		FileSearcher.main(new String[] { "foldername=" + input, "writefn=" + fastQFiles,
				"searchstrings=" + fastQSearchStrings, "forbiddenstrings=" + forbiddenSearchStrings });
		return fastQFiles;
	}
}
