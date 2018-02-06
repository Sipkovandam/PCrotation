package GeneNetwork;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;

import Tools.FileUtils;
import Tools.Script;

public class ReactomeTermRetriever extends Script<ReactomeTermRetriever>
{
	//uses the gene-network api to get 1 reactome term per gene (based on co-expressed genes enrichment)
	String geneNamesFn = null;//must contain all gene names to be queried in the first column
	String url = "http://molgenis27.target.rug.nl/api/v1/gene/";
	
	@Override
	public void run()
	{
		try
		{
			String writeFn = FileUtils.removeExtention(geneNamesFn)+"_ReactomeTerms.txt";
			
			BufferedReader geneReader = FileUtils.createReader(geneNamesFn);
			BufferedWriter reactomeWriter = FileUtils.createWriter(writeFn);
			
			//reads the first line and writes it + the extra columns
			writeHeader(geneReader, reactomeWriter);
			
			String line = null;
			int l = 0;
			while((line = geneReader.readLine())!=null)
			{
				String gene = line.split("\t")[0];
				String webRequest = url+gene;
				System.out.println(webRequest);
				String reactomeTerm = getUrlSource(webRequest, gene);
		
				//System.out.println(reactomeTerm);
				
				if(reactomeTerm!=null)
					reactomeWriter.write(line + "\t" + reactomeTerm + "\n");
				else
					reactomeWriter.write(line + "\t" + "" +"\t"+ "" + "\n");
				l++;
			}
			reactomeWriter.close();
		}catch(Exception e){e.printStackTrace();}
	}
	
	private void writeHeader(BufferedReader geneReader, BufferedWriter reactomeWriter) throws IOException
	{
		String header = geneReader.readLine();
		reactomeWriter.write(header+"\tReactome term\tp-Value\n");
	}

	private static String getUrlSource(String url, String gene) throws IOException {
        URL geneNetwork = new URL(url);
        URLConnection yc =geneNetwork.openConnection();
    	BufferedReader in = null;
        try
        {
        	in = new BufferedReader(new InputStreamReader(
            yc.getInputStream(), "UTF-8"));
        }catch(Exception e)
        {
        	System.out.println("Gene not present in database:\t" + gene);
        	return null;
        }
        
        String webLine;
        StringBuilder a = new StringBuilder();
        while ((webLine = in.readLine()) != null)
        {
        	if(webLine.contains("REACTOME:"))
        	{
        		String firstSplit = webLine.split("REACTOME:")[1];
        		String term = firstSplit.split("\"")[0];
        		String zScore = in.readLine().split("\"zScore\": ")[1].split(",")[0];
        		String pValue = in.readLine().split("\"pValue\": ")[1].split(",")[0];
        		in.close();
        		return term+"\t"+pValue;
        	}
            //a.append(webLine);
        }
        in.close();

        return null;//a.toString();
    }
	
	
}
