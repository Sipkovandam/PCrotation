package PrepData;

import java.io.IOException;

import Analyses.WilcoxonMannWhitney;
import MatrixScripts.MatrixStruct;
import umcg.genetica.math.stats.Correlation;

public class SampleCorrelation {

	public static void main(String[] args) throws IOException 
	{
		String samplesFN = "E:/Groningen/Data/Monogenetic_disease_samples_RadBoud/directPCA_Rlog_0.2/18DownSyndrome26Normal2Cancer_counts/centered.txt";
		String publicFN= "E:/Groningen/Data/Monogenetic_disease_samples_RadBoud/directPCA_Rlog_0.2/MATRIX_Centered.txt";
		String writeFN = "E:/Groningen/Data/Monogenetic_disease_samples_RadBoud/directPCA_Rlog_0.2/18DownSyndrome26Normal2Cancer_counts/correlationWithPublic.txt";
		
		if(args.length< 1)
			checkArgs(args);
		
		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase()){
				case "filename":
					samplesFN =value;
					break;
				case "publicfn":
					publicFN = value;
					break;
				case "writefn":
					writeFN = value;
					break;
				default:
					checkArgs(args);
					System.out.println("Incorrect argument supplied; exiting");
					System.exit(1);
			}
		}
		
		getCorrelation(samplesFN,publicFN, writeFN);		
	}


	private static void getCorrelation(String samplesFN, String publicFN, String writeFN) throws IOException 
	{
		MatrixStruct samples = new MatrixStruct(samplesFN);
		MatrixStruct publicSamples = new MatrixStruct(publicFN);
		samples.putGenesOnRows();
		publicSamples.putGenesOnRows();
		samples.keepRows(publicSamples);
		
		MatrixStruct correlations = new MatrixStruct(publicSamples.cols(),samples.cols());
		correlations.setColHeaders(samples.getColHeaders());
		correlations.setRowHeaders(publicSamples.getColHeaders());

		for(int c = 0; c < samples.cols(); c++)
		{
			System.out.println("Sample " + c + "/" + samples.cols());
			double[] sampleVals = samples.getColValues(c);
			for(int cP = 0; cP < publicSamples.cols(); cP++)
			{
				//System.out.println("PublicSample" + cP + "/" + publicSamples.cols());
				double[] publicSampleVals = publicSamples.getColValues(cP);
				double correlation = Correlation.correlate(sampleVals, publicSampleVals);
				correlations.matrix.set(cP,c,correlation);
			}
		}
		correlations.write(writeFN);
		System.out.println("File written to: " + writeFN);
	}
	static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
			return;
		System.out.println("Wrong arguments");
		System.exit(1);
	}
}