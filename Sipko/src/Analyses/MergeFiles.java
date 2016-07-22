package Analyses;

import java.io.IOException;

import pca.MatrixStruct;

public class MergeFiles {

	public static void main(String[] args) throws IOException 
	{
		boolean merge = false; //if false it only sorts the fn1 matrix according to the fn2 matrix
		String fn1="E:/Groningen/Data/PublicSamples/05-2016/Tessa/est_counts_nocancernocellline_Rlog_covariance_Chr1-22_add0.0BeforeGeoMean/18DownSyndrome26Normal2Cancer_counts/1healhtySampleCompare.txt";
		String fn2="E:/Groningen/Data/PublicSamples/05-2016/Tessa/est_counts_nocancernocellline_Rlog_covariance_Chr1-22_add0.0BeforeGeoMean/18DownSyndrome26Normal2Cancer_counts/1DownSampleCompare.txt";

//		String fn1="E:/Groningen/Data/PublicSamples/05-2016/Tessa/est_counts_nocancernocellline_Rlog_covariance_Chr1-22_add0.0BeforeGeoMean/CountsGENES_5GPM_46Samples/centered_transposed.txt";
//		String fn2="E:/Groningen/Data/PublicSamples/05-2016/Tessa/est_counts_nocancernocellline_Rlog_covariance_Chr1-22_add0.0BeforeGeoMean/CountsGENES_5GPM_46Samples/PC_1-300_DevidedBySTdevs.txt";

		
		MatrixStruct f1= new MatrixStruct(fn1);
		MatrixStruct f2= new MatrixStruct(fn2);

		if(merge)
			f1=f1.mergeColumns(f2);
		else
			f2.keepRows(f1);
		
		f1.write(fn1.replace(".txt", "").replace(".gz", "")+"_merged.txt");
	}

}
