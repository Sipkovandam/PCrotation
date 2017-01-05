package STAR;

public class SpliceCount 
{

	private String sampleName;
	private boolean includeInReference;
	private int readsOverlapping;
	private double[] relativeAbundances;
	
	SpliceCount(String sampleName, int readsOverlapping, double[] relativeAbundances, boolean includeInReference)
	{
		this.sampleName = sampleName;
		this.readsOverlapping = readsOverlapping;
		this.relativeAbundances = relativeAbundances;
		this.includeInReference = includeInReference;
	}
	
	/**
	 * @return the includeInReference
	 */
	public boolean getIncludeInReference() {
		return includeInReference;
	}
	
	
	/**
	 * @return the sampleName
	 */
	public String getSampleName() {
		return sampleName;
	}
	/**
	 * @param sampleName the sampleName to set
	 */
	public void setSampleName(String sampleName) {
		this.sampleName = sampleName;
	}
	/**
	 * @return the relativeAbundances
	 */
	public double[] getRelativeAbundances() {
		return relativeAbundances;
	}
	/**
	 * @param relativeAbundance the relativeAbundance to set
	 */
	public void setRelativeAbundance(double[] relativeAbundances) {
		this.relativeAbundances = relativeAbundances;
	}
	/**
	 * @return the readsOverlapping
	 */
	public int getReadsOverlapping() {
		return readsOverlapping;
	}
	/**
	 * @param readsOverlapping the readsOverlapping to set
	 */
	public void setReadsOverlapping(int readsOverlapping) {
		this.readsOverlapping = readsOverlapping;
	}
	/**
	 * @return the includeInReference
	 */
	public boolean isIncludeInReference() {
		return includeInReference;
	}
	/**
	 * @param includeInReference the includeInReference to set
	 */
	public void setIncludeInReference(boolean includeInReference) {
		this.includeInReference = includeInReference;
	}
	
	
}
