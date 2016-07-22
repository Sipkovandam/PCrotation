package PCA;

import java.io.IOException;

public class zTransform {

	static public String FN = "";
	static public String writeFN = null;
	static public boolean rows = true;
	
	public static void main(String[] args) throws IOException 
	{
		
		if(args.length == 0)
			checkArgs(args);
		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase()){
				case "filename":
					FN =value;
					break;
				case "writefn":
					writeFN = value;
					break;
				case "rows":
					rows = Boolean.parseBoolean(value);
					break;
				default:
					checkArgs(args);
					System.out.println("Incorrect argument supplied; exiting");
					System.exit(1);
			}
		}
		if(writeFN == null)
		{
			if(rows)
				writeFN = FN.replace(".txt", ".zTransformed_Rows.txt");
			if(!rows)
				writeFN = FN.replace(".txt", ".zTransformed_Cols.txt");
		}

		MatrixStruct sample = new MatrixStruct(FN);
				MatrixStruct stDevs = null;
		if(rows)
			stDevs=sample.stDevRows();
		else
			stDevs=sample.stDevCols();
		
		pca.PCA.log("Writing STdevs " + writeFN.replace(".txt", ".Stdevs.txt"));
		stDevs.write(writeFN.replace(".txt", ".Stdevs.txt"));
		pca.PCA.log("Divide all gene values by STdev for each sample");	
		sample.divideBy(stDevs,rows);//false corrects columns, true corrects rows
		
		pca.PCA.log("Writing zTransformed matrix" + writeFN);
		sample.write(writeFN);
		pca.PCA.log("Done");
	}
	static void checkArgs(String[] args) 
	{
		System.out.println("Wrong arguments, needs at least 1 argument:\n"
				+ "fileName=<filename.txt> - Filename of the file you want zTransform (set stdevs to 1)\n"
				+ "writeFN=<writeFN.txt> - Write filename (default=<filename.txt>.replace(.txt, .ztransformed_<rows/cols>.txt)\n"
				+ "rows=<true/false> - if true sets stdev of each row to 1, otherwise cols (default=true)");
		System.exit(1);
	}
}
