package STAR;

import java.util.HashMap;

public class SpliceSitesCount 
{
	//private String spliceSite = null;
	private int maxOverhang;//maximum size of the overhang of this junction detected in all samples
	private int reads;//total number of reads overlapping this junction in all samples together
	private int samples;//number of samples this junction is observed in
	private int samplesAboveCutoff;//number of reference samples this junction is observed in with 8 (2^3) or more reads overlapping this junction
	private String annotated;
	
	SpliceSitesCount(String annotated)
	{
		this.setAnnotated(annotated);
	}
	
	SpliceSitesCount(int overhang, int reads, int samplesAboveCutoff, String annotated)
	{
		this.setMaxOverhang(overhang);
		this.setReads(reads);
		this.samples =1;
		this.samplesAboveCutoff=samplesAboveCutoff;
		this.setAnnotated(annotated);
	}

	public void incrementSamples()
	{
		samples++;
	}
	public int getSamples()
	{
		return samples;
	}
	
	public void incrementSamplesAboveCutoff()
	{
		samplesAboveCutoff++;
	}
	public int getSamplesAboveCutoff()
	{
		return samplesAboveCutoff;
	}
	
	public int getMaxOverhang() {
		return maxOverhang;
	}

	public void setMaxOverhang(int maxOverhang) {
		this.maxOverhang = maxOverhang;
	}

	public String getAnnotated() {
		return annotated;
	}

	public void setAnnotated(String annotated) {
		this.annotated = annotated;
	}

	public int getReads() {
		return reads;
	}

	public void setReads(int reads) {
		this.reads = reads;
	}
	
	public void incrementReads(int reads) {
		this.setReads(this.getReads()+reads);
	}
}
