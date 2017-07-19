package PCA;

import java.io.IOException;
import java.util.ArrayList;

import MatrixScripts.MyMatrix;
import umcg.genetica.math.stats.Correlation;

public class RemoveDuplicates {
	//Removes duplicate samples. Duplicates are defined as <removeDuplicates> correlated.
	//Warning: Assumes duplicates are next to eachother (this was faster and is usually the case in biological data, but definetly not always)
	
	public static void main(String[] args) throws IOException 
	{
		String expressionFN = "";
		String writeFolder = "";
		boolean writeAll = true;
		double removeDuplicates = 1;	

		for(int a = 0; a < args.length; a++)
		{
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase()){
				case "filename":
					expressionFN =value;
					break;
				case "writefolder":
					writeFolder = value;
					break;
				case "writeall":
					writeAll = Boolean.parseBoolean(value);
					break;
				case "removeduplicates":
					removeDuplicates = Double.parseDouble(value);
					break;
				default:
					checkArgs(args);
					System.out.println("Incorrect argument supplied; exiting");
					System.exit(1);
			}
		}
		
		MyMatrix expressionStruct = new MyMatrix(expressionFN);
		removeDuplicates(expressionStruct, removeDuplicates, writeFolder,writeAll);
	}

	public static void removeDuplicates(MyMatrix expression, double duplicateCutoff, String writeFolder, boolean write) throws IOException {
		//this function assumes duplicate rows are always next to each other (as is usually the case)
		//This saves some computational time
		if(duplicateCutoff >= 1)
			return;
		JuhaPCA.PCA.log(" 2. Removing duplicates (r>"+duplicateCutoff+")");
		ArrayList<Integer> rowsToRemove = new ArrayList<Integer>();
		for(int r = 0; r < expression.rows()-1; r++)
		{
			double correlation = Correlation.correlate(expression.getRowValues(r), expression.getRowValues(r+1));
			if(correlation > duplicateCutoff)
				rowsToRemove.add(r+1);
		}
		System.out.println("Removing  " + rowsToRemove.size() + " duplicates");
		MyMatrix adjustedMatrix = new MyMatrix(expression.rows()-rowsToRemove.size(),expression.cols());
		adjustedMatrix.setColHeaders(expression.getColHeaders());
		int a = 0;
		int outRow = 0;
		if(rowsToRemove.size()==0)
			return;
		for(int r = 0; r < expression.rows(); r++)
		{
			int skip = rowsToRemove.get(a);
			if(skip == r)
			{
				if(a<rowsToRemove.size()-1)
					a++;
				continue;
			}
			adjustedMatrix.setRow(outRow, expression.getRowHeaders()[r], expression.getRowValues(r));
			outRow++;
		}
		String duplicatesRemovedFN = writeFolder+"DuplicatesRemoved.txt";
		if(write)adjustedMatrix.write(duplicatesRemovedFN);
		expression = adjustedMatrix;
	}
	static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
			return;
		System.out.println("Wrong arguments");
		System.exit(1);
	}
}
