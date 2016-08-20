package PCA;

public class LowestToAverage {
	//sets the lowest value in a dataset to the average
	
	public static void main(String[] args) 
	{
		String expressionFN = "";
		String writeFolder = "";
		boolean writeAll = true;
		
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
				default:
					checkArgs(args);
					System.out.println("Incorrect argument supplied; exiting");
					System.exit(1);
			}
		}
		
		MatrixStruct expressionStruct = new MatrixStruct(expressionFN);
		lowestToAverage(expressionStruct);
	}

	public static void lowestToAverage(MatrixStruct expressionStruct) 
	{
		for(int r = 0; r < expressionStruct.rows(); r++)
		{
			double lowest = Double.MAX_VALUE;
			double sum = 0;
			for(int c = 0; c < expressionStruct.cols(); c++)
			{
				double val = expressionStruct.matrix.get(r, c);
				sum+= val;
				if(val<lowest)
					lowest = val;
			}
			double average = sum/expressionStruct.cols();
			for(int c = 0; c < expressionStruct.cols(); c++)
			{
				if(expressionStruct.matrix.get(r, c) == lowest)
					expressionStruct.matrix.set(r, c, average);
			}
		}
	}
	static void checkArgs(String[] args) 
	{
		if(System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
			return;
		System.out.println("Wrong arguments");
		System.exit(1);
	}

}
