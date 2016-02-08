package PCA;

public class Transpose 
{

	public static void main (String[] args)
	{
		String fileName = "E:/Groningen/Data/PublicSamples/Test8/est_counts_nocancernocelllineSTdevRND/GENE.eigenvectorsSmall.txt";
			
		Matrix matrix = new Matrix(fileName);
		matrix.transpose();
		matrix.write(fileName.replace(".txt", "_transposed.txt"));
//		matrix.write(fileName.replace(".txt", "_first10geneEVs.txt"), -1,10);
//		matrix.write(fileName.replace(".txt", "_first100geneEvs.txt"), -1,100);
//		matrix.write(fileName.replace(".txt", "_transposed_100.txt"), 100, -1);
//		matrix.write(fileName.replace(".txt", "_transposed_300.txt"), 300, -1);
//		matrix.write(fileName.replace(".txt", "_transposed_500.txt"), 500, -1);
//		matrix.write(fileName.replace(".txt", "_transposed_1000.txt"), 1000, -1);
		System.out.println("Done");
	}
}
