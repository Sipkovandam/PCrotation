package PrepData;

import java.util.ArrayList;
import java.net.*;
import java.io.*;

public class GetStudies 
{
	public static void main(String[] args) throws IOException
	{
		String[] studyNames = new String[]{"PRJNA227344","PRJNA213398","SRP011927"};
		String writeFN = "E:/Groningen/Data/PublicSamples/Remapping samples/samplesToGet.txt";
		ArrayList<String> webLinks = new ArrayList<String>();
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(writeFN)));
		for(int s = 0; s < studyNames.length; s++)
		{
			//study_accession	secondary_study_accession	sample_accession	secondary_sample_accession	experiment_accession	run_accession	submission_accession	tax_id	scientific_name	instrument_platform	instrument_model	library_name	library_layout	nominal_length	library_strategy	library_source	library_selection	read_count	base_count	center_name	first_public	experiment_title	study_title	study_alias	experiment_alias	run_alias	fastq_bytes	fastq_md5	fastq_ftp	fastq_aspera	fastq_galaxy	submitted_bytes	submitted_md5	submitted_ftp	submitted_aspera	submitted_format	submitted_galaxy	sra_bytes	sra_md5	sra_ftp	sra_aspera	sra_galaxy	cram_index_ftp	cram_index_aspera	cram_index_galaxy	col_tax_id	col_scientific_name	sample_alias	broker_name
			//http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=PRJEB1195&result=read_run&fields=study_accession%2Cfastq_ftp%2Cfastq_md5%2Cfastq_bytes%2Ctax_id

			String webLink = "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession="+ studyNames[s]+"&result=read_run&fields=study_accession%2Csecondary_study_accession%2Csample_accession%2Csecondary_sample_accession%2Cexperiment_accession%2Crun_accession%2Csubmission_accession%2Ctax_id%2Cscientific_name%2Cinstrument_platform%2Cinstrument_model%2Clibrary_name%2Clibrary_layout%2Cnominal_length%2Clibrary_strategy%2Clibrary_source%2Clibrary_selection%2Cread_count%2Cbase_count%2Ccenter_name%2Cfirst_public%2Cexperiment_title%2Cstudy_title%2Cstudy_alias%2Cexperiment_alias%2Crun_alias%2Cfastq_bytes%2Cfastq_md5%2Cfastq_ftp%2Cfastq_aspera%2Cfastq_galaxy%2Csubmitted_bytes%2Csubmitted_md5%2Csubmitted_ftp%2Csubmitted_aspera%2Csubmitted_format%2Csubmitted_galaxy%2Csra_bytes%2Csra_md5%2Csra_ftp%2Csra_aspera%2Csra_galaxy%2Ccram_index_ftp%2Ccram_index_aspera%2Ccram_index_galaxy%2Ccol_tax_id%2Ccol_scientific_name%2Csample_alias%2Cbroker_name";
			System.out.println(webLink);
			getWebPage(webLink, writer,s);
		}
		writer.close();
	}
	public static void getWebPage(String webLink, BufferedWriter writer, int s) throws IOException
	{
		URL url = new URL(webLink);
	
		BufferedReader in = new BufferedReader(new InputStreamReader(url.openStream()));
		String inputLine;
		if(s>0)
			in.readLine();
		while ((inputLine = in.readLine()) != null) {
		  writer.write(inputLine + "\n");
		}
		
		in.close();
	}
}
