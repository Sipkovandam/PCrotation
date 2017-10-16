package TextEditing;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import Tools.FileUtils;
import Tools.Script;

public class FileSplitter extends Script<FileSplitter>
{
	//splits a file with the gene names in a column based on gene names into multiple files
	//optionally the names in the column can be further split by e.g. "__"
	//this script uses a large amount of memory; Approximately 240 GB for 1 million splice sites (Array Lists and Hashtables aren't very memory efficient...)
	String fnComment = "File to be split based on genes. Has gene names in first column";
	String fn = null;
	String spliceNamesColComment = "Column in which the splice file names are located (lines in new file include this and all columns behind it). This does not start from [0] index, so 1 is the 1st/rowname column";
	int spliceNamesCol = 2;
	String geneNameColumnComment = "Column in which the gene names are located (separated by <splitRegex>)";
	int geneNameColumn = 0;
	String splitFirstColumnInsteadComment = "If true, splits the first column based on splitRegex instead, to obtain the gene name";
	boolean splitFirstColumnInstead = false;

	String splitRegexComment = "splits the rowname based on this and takes the [0] index as filename to write to; retains the whole line;";
	String splitRegex = ",";//splits the rowname based on this and takes the [0] index as filename to write to; retains the whole line;
	String writeFolderComment = "Folder to write the results to";
	String writeFolder = null;
	String gzipFilesComment = "If true, gzips the output files";
	boolean gzipFiles = true;

	public void run()
	{
		try
		{
			if(writeFolder==null)
				writeFolder=new File(fn).getParent()+"/perGene/";
			if(new File(writeFolder).exists() && new File(writeFolder).list().length>1)
			{
				p("WriteFolder needs to be empty:\nrm -r " + writeFolder);
				p("Exiting");
				System.exit(2);
			}
				
			FileUtils.makeDir(writeFolder);

			BufferedReader reader = FileUtils.createReader(fn);
			String header = null;
			
			if(splitFirstColumnInstead)
				header=reader.readLine();
			else
				header=reader.readLine().split("\t",spliceNamesCol+1)[spliceNamesCol];
			
			String line = null;
			BufferedWriter splicePerGeneWriter = null;
			String currentGene = "";
			String currentGeneFn = "";
			int n = 0;

			String writeExtencion = ".txt";
			if(gzipFiles)
				writeExtencion = ".txt.gz";
			
			HashMap<String, BufferedWriter> writers=new HashMap<String, BufferedWriter>();
			ArrayList<String> geneWriterOrder = new ArrayList<>();
			
			//put all the genes and corresponding writeLines into a hash.
			HashMap<String, ArrayList<String>> geneToExpressionLines = new HashMap<String, ArrayList<String>>();
			while ((line = reader.readLine()) != null)
			{

				if (n % 100000 == 0)
					p(n + " lines read");

				String[] genes = null;
				String[] eles = line.split("\t");
				if (splitFirstColumnInstead)
				{
					genes = new String[] { eles[0].split(splitRegex)[0] };//splits the first column by <splitRegex> and gets the first element to get the gene name
				}
				else
					genes = eles[geneNameColumn].split(splitRegex);

				//get all the spliceLines that should be written for this gene
				for (String gene : genes)
				{
					BufferedWriter geneWriter = getWriter(gene, writers, geneWriterOrder, writeExtencion,header);
					
//					ArrayList<String> geneLines = geneToExpressionLines.get(gene);
//					if (geneLines == null)
//						geneLines = new ArrayList<String>();

					StringBuilder lineToAdd = new StringBuilder();
					if(splitFirstColumnInstead)
						lineToAdd.append(line);
					else
						lineToAdd.append(gene).append("__").append(line.split("\t", spliceNamesCol+1)[spliceNamesCol]);
					lineToAdd.append("\n");
					geneWriter.write(lineToAdd.toString());
//					geneLines.add(lineToAdd);
//					
//					geneToExpressionLines.put(	gene,
//												geneLines);
				}
				n++;
			}
			//close all remaining writers
			closeWriters(writers);

//it doesnt fit in memory.
			//write the writeLines per gene
//			String finalHeader = header;
//			geneToExpressionLines.forEach((	gene,
//											writeLines) ->
//			{
//				try
//				{
//					String writeExtencion = ".txt";
//					if(gzipFiles)
//						writeExtencion = ".txt.gz";
//					BufferedWriter geneWriter = FileUtils.createWriter(writeFolder + gene + writeExtencion);
//					geneWriter.write(finalHeader + "\n");
//					for (String writeLine : writeLines)
//					{
//						geneWriter.write(writeLine + "\n");
//					}
//					geneWriter.close();
//				} catch (Exception e)
//				{
//					e.printStackTrace();
//				}
//			});

			p("Finished! Files written at:\t" + writeFolder);
		} catch (Exception e)
		{
			e.printStackTrace();
		}

	}

	private void closeWriters(HashMap<String, BufferedWriter> writers)
	{
		writers.forEach((g,v)->{try{v.close();} catch (Exception e){}});
	}

	//keeps a miximum of 100 writers open at the same time
	private BufferedWriter getWriter(	String gene,
										HashMap<String, BufferedWriter> writers,
										ArrayList<String> geneWriterOrder,
										String writeExtencion, String header) throws FileNotFoundException, IOException
	{
		BufferedWriter writer = writers.get(gene);
		String fileName = writeFolder + gene + writeExtencion;
		int maxWritersOpen =5;
		if(writer==null)
		{
			if(new File(fileName).exists())
				writer=FileUtils.createWriter(fileName,true);
			else
			{
				writer=FileUtils.createWriter(fileName);
				writer.write(header+"\n");
			}
				
			writers.put(gene, writer);
			geneWriterOrder.add(gene);
			if(geneWriterOrder.size()>maxWritersOpen)
			{
				String oldGene=geneWriterOrder.get(0);
				writers.get(oldGene).close();
				writers.remove(oldGene);
				geneWriterOrder.remove(0);
			}
		}

		return writer;
	}

	public String getFn()
	{
		return fn;
	}

	public void setFn(String fn)
	{
		this.fn = fn;
	}

	public int getSpliceNamesCol()
	{
		return spliceNamesCol;
	}

	public void setSpliceNamesCol(int spliceNamesCol)
	{
		this.spliceNamesCol = spliceNamesCol;
	}

	public String getSplitRegex()
	{
		return splitRegex;
	}

	public void setSplitRegex(String splitRegex)
	{
		this.splitRegex = splitRegex;
	}

	public String getWriteFolder()
	{
		return writeFolder;
	}

	public void setWriteFolder(String writeFolder)
	{
		this.writeFolder = writeFolder;
	}

	public boolean isGzipFiles()
	{
		return gzipFiles;
	}

	public void setGzipFiles(boolean gzipFiles)
	{
		this.gzipFiles = gzipFiles;
	}

	public int getGeneNameColumn()
	{
		return geneNameColumn;
	}

	public void setGeneNameColumn(int geneNameColumn)
	{
		this.geneNameColumn = geneNameColumn;
	}

	public boolean isSplitFirstColumnInstead()
	{
		return splitFirstColumnInstead;
	}

	public void setSplitFirstColumnInstead(boolean splitFirstColumnInstead)
	{
		this.splitFirstColumnInstead = splitFirstColumnInstead;
	}

	public String getGzipFilesComment()
	{
		return gzipFilesComment;
	}

	public void setGzipFilesComment(String gzipFilesComment)
	{
		this.gzipFilesComment = gzipFilesComment;
	}
}
