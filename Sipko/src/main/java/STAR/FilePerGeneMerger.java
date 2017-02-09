package STAR;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import PCA.MatrixString;
import Tools.FileSearcher;
import Tools.FileUtils;
import Tools.Script;
import umcg.genetica.containers.Pair;

public class FilePerGeneMerger extends Script<FilePerGeneMerger> 
{
	//also transposes
	//rownames need to be in same order (so before transposin)
	//mergerges per gene output files of spliceMerger and ExonexpressionMerger
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 8465877177124291824L;
	public String folderNameComment = "/home/directory/; Foldername in which the files are located that should be merged";
	public String folderName = "E:/Groningen/Splicing/100BPcap_analysis/Results/PerGene/Ratios/";
	public String writeFnComment = "/home/directory/result.txt; Filename of the merged file (result)";
	public String writeFn = "E:/Groningen/Splicing/100BPcap_analysis/Results/PerGene/ratiosMerged.txt";
	public transient ArrayList<String> filesToMerge = null;
	
	public FilePerGeneMerger()
	{
	}
	
	public FilePerGeneMerger(ArrayList<String> filesToMerge, String writeFn)
	{
		this.filesToMerge = filesToMerge;
		this.writeFn=writeFn;
	}
	
	@Override
	public void run()
	{
		try
		{
			if(filesToMerge==null)
			{
				FileSearcher fileSearcher = new FileSearcher(new File(this.folderName).getParent());
				fileSearcher.setFolders(this.folderName);
				fileSearcher.setSearchStrings(new String[]{".txt"});
				String fnsFn = new File(this.folderName).getParent()+"/Fns.txt";
				fileSearcher.setWriteName(fnsFn);
				fileSearcher.run();
				filesToMerge=FileUtils.readArrayList(fnsFn);
			}
			if(writeFn==null)
			{
				p("writeFn variable not set; exiting");
				System.exit(2);
			}
			
			Set<String> writtenSplices= new HashSet<>();
			BufferedWriter mergedWriter = FileUtils.createWriter(writeFn);
			String[] headers = null;
			for(String fn : filesToMerge)
			{
				MatrixString matrix = new MatrixString(fn);
				p(fn);
				matrix.transpose();
				List<Object> retVals = checkColNames(headers, matrix, mergedWriter);
				boolean sameColNameOrder = (boolean) retVals.get(0);
				headers = (String[]) retVals.get(1);
				for(int r = 0; r < matrix.rows(); r++)
				{
					if(sameColNameOrder)
					{
						writtenSplices.add(matrix.rowNames[r]);
						mergedWriter.write(matrix.getAsLine(matrix.values[r], matrix.rowNames[r]));
					}
					else
					{
						p("Colnames do not have same order. No function was implemented to manage this yet; exiting");
						System.exit(2);
					}
				}
			}
			mergedWriter.close();
			p("Done merging files. File written to: " + writeFn);
		}catch(Exception e){e.printStackTrace();}
	}

	private List<Object> checkColNames(String[] headers, MatrixString matrix, BufferedWriter mergedWriter) throws IOException {
		boolean returnVal = true;
		if(headers == null)
		{
			headers= matrix.colNames;
			mergedWriter.write(matrix.getAsLine(headers,""));
		}
		else
		{
			for(int c = 0; c < matrix.colNames.length; c++)
			{
				if(!matrix.colNames[c].equals(headers[c]))
				{

					returnVal=false;
				}
			}
		}
		return r(returnVal,headers);
	}
}
