package STAR;

import java.util.ArrayList;

import Slurm.ClusterHandler;
import Tools.FileSearcher;
import Tools.FileUtils;
import Tools.Script;

public class FastqDeduplicator extends Script<FastqDeduplicator>
{	
	FileSearcher filesToDeduplicateSearcher = new FileSearcher();
	ClusterHandler<FastqDeduplicator_ClusterHandler> slurm = new ClusterHandler<FastqDeduplicator_ClusterHandler>("BBMap/35.69-Java-1.7.0_80");
	
	private String pairStringsComment = "[\"_R1_\",\"_R2_\",\"_1.fq\",\"_2.fq\"]; MANDATORY for paired end data // comma separated list of string pairs defining the difference between forward read and backward read files.  For example, immagine a file name fastqFile_1.fq - to obtain the complementary file in this file name _R1_ is replaced iwth _R2_ and _1.fq is replaced with _2.fq. Since _R1_ is not present in the file name but _1.fq is, the complementary file becomes fastqFile_2.fq and these 2 files then are used by STAR. If you files do not contain any of these strings STAR will map the data as if it was single end data";
	private String[] pairStrings = new String[] { "_R1", "_R2", "_1.fq", "_2.fq" };
	
	@Override
	public void run()
	{
		try
		{
			filesToDeduplicateSearcher.run();
			slurm.setFastQFiles(filesToDeduplicateSearcher.getWriteName());
			slurm.run();
			
		}catch(Exception e){e.printStackTrace();}
		
	}

	public String[] getPairStrings()
	{
		return pairStrings;
	}

	public void setPairStrings(String[] pairStrings)
	{
		this.pairStrings = pairStrings;
	}

	public FileSearcher getFilesToDeduplicateSearcher()
	{
		return filesToDeduplicateSearcher;
	}

	public void setFilesToDeduplicateSearcher(FileSearcher filesToDeduplicateSearcher)
	{
		this.filesToDeduplicateSearcher = filesToDeduplicateSearcher;
	}
}
