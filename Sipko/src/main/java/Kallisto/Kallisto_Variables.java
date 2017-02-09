package Kallisto;

import java.io.Serializable;

public class Kallisto_Variables implements Cloneable, Serializable
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	private String kallistoIndexFileComment = "/root/directory/hg19.v75.cdna.all.42.2.idx; MANDATORY //Annotation file Kallisto should use. This file has to be created by Kallisto prior to running this script";
	private String kallistoIndexFile = "/groups/umcg-wijmenga/tmp04/umcg-svandam/Data/RNAseq/Annotation/hg19.v75.cdna.all.42.2.idx";

	public String getKallistoIndexFile() {
		return kallistoIndexFile;
	}

	public void setKallistoIndexFile(String kallistoIndexFile) {
		this.kallistoIndexFile = kallistoIndexFile;
	}

}
