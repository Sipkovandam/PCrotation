package MatrixScripts;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.LineNumberReader;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Date;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.stream.Stream;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.nio.file.Path;

public class MatrixString
{
	//a matrix class that is not limited to max java array size
	
	public String firstField;
	public String[] colNames;
	public String[] rowNames;
	public String[][] values;
	public boolean verbose = false;
	private Hashtable<String,Integer> rowHash = null;
	private Hashtable<String,Integer> colHash = null;
	protected static final String ENCODING = "ISO-8859-1";
	public static final int DEFAULT_BUFFER_SIZE = 4096;
	public class GetVal
	{
		String[][] values;
		
		GetVal(String[][] matrix)
		{
			values = matrix;
		}
		public String get(int r, int c)
		{
			return values[r][c];
		}
		public void set(int r, int c, String d) {
			values[r][c]=d;
		};
	}
	public GetVal matrix = new GetVal(values);
	
	public MatrixString()
	{
		
	}
	public MatrixString(int x, int y)
	{
		rowNames = new String[x];
		colNames = new String[y];
		values = new String[x][y];
		matrix = new GetVal(values);
		
		numberNames(rowNames, "Row");
		numberNames(colNames, "Col");		
	}
	
	public MatrixString(int x, int y, String[] colNames, String[] rowNames)
	{
		this.rowNames = rowNames;
		this.colNames = colNames;
		values = new String[x][y];
		matrix = new GetVal(values);
	}
	
	
	public MatrixString(String fileName)
	{
		readFile(fileName,true,true);
	}
	public MatrixString(Path fileName)
	{
		readFile(fileName.toString(),true,true);
	}
	public MatrixString(String fileName,int x, int y)//read file with a maximum of xRows and yCols
	{
		readFile(fileName,true,true, x, y);
	}
	
	public MatrixString(String fileName, boolean hasRowColNames)
	{
		readFile(fileName,false,false);
	}
	public MatrixString(String[] rowNames, String[] colNames, String[][] matrix) 
	{
		this.rowNames = rowNames;
		this.colNames = colNames;
		this.values = matrix;
		this.matrix = new GetVal(values);
	}

	public void readFile(String fileName)
	{
		readFile(fileName,true,true);
	}
	public void readFile(String fileName,boolean hasRowNames, boolean hasColName)
	{
		readFile(fileName, hasRowNames, hasColName, -1, -1);
	}
	
	private void numberNames(String[] names, String extra) 
	{
		for(int x = 0; x < names.length; x++)
		{
			names[x] = extra+ "" + Integer.toString(x);
		}
	}
	
	public void readFile(String fileName,boolean hasRowNames, boolean hasColNames,int maxX, int maxY)
	{
		int nRows = 0;
		if(maxX<1)
			nRows = getRowNumber(fileName)-1;
		else
			nRows = maxX;
				
		try 
		{
			BufferedReader reader = getReader(fileName);
			
			String line = null;
			line = reader.readLine();
			String[] eles = line.split("\t");
			int nCols= eles.length;
			if(nCols > maxY && maxY > 0)
				nCols = maxY;
			
			if(verbose)
			{
				System.out.println("Input file has the following format (rows,columns): (" + nRows + "," + nCols + ")");
				//System.out.println("This includes 1 row and 1 column for the row/col names");
			}	
			
			//Matrix looks like this:
			//X n n n
			//n v v v
			//n v v v
			//n v v v
			//Where "n" is row or column name and "v" is value
			
			String secondLine = reader.readLine(); 
			int secondLineCols= secondLine.split("\t").length;
			int minusCol = 1;
			if(secondLineCols == nCols+1)//first line is missing a cell
				minusCol = 0;
			
			//Put the colNames into the matrix
			firstField = eles[0];
			if(hasColNames && hasRowNames)
				values = new String[nRows-1][nCols-minusCol];
			else if (hasColNames && !hasRowNames)
				values = new String[nRows][nCols-minusCol];
			else if (!hasColNames && hasRowNames)
				values = new String[nRows-1][nCols];
			else
				values = new String[nRows][nCols];
			matrix = new GetVal(values);
			if(hasColNames)
			{
				colNames = new String[nCols-minusCol];//-1 because one of them is for the rowNames, that column does not have a column name
				//Deal with first line missing first cell
				for(int y = 0; y < nCols-minusCol; y++)
				{
					colNames[y] = eles[y+minusCol];
				}
			}
			else //reset reader to start of file, probably not the best way of doing it...
			{
				colNames = new String[nCols];//-1 because one of them is for the rowNames, that column does not have a column name
				
				reader.close();
				for(int y = 0; y < eles.length; y++)
				{
					colNames[y] = "Col"+Integer.toString(y);
				}
				reader = getReader(fileName);
			}
			if(hasRowNames)
				rowNames = new String[nRows-1];
			else
				rowNames = new String[nRows];
		
			int num = 0;//number to substract if it has rowNames (1) or if it does not (0)
			if(hasRowNames)
				num = 1;
			
			int x = 0;
			if(hasColNames == true)
				line = secondLine;
			else
				line = reader.readLine();
			
			//System.out.println("");
			double nX = 0;
			double blockSize = 1000000000;//print after 100.000.000 cells have been read
			int printPoint = (int)(blockSize/nCols);
			
			DecimalFormat f = new DecimalFormat("#.##");
			do
			{
				//System.out.println(line);
				if(x >= rowNames.length)// in case the last line of the file is corrupt (e.g. an incomplete tranfer) the last line will not be included
				{
					break;
				}
				
				if(nX % printPoint == 0 && nX != 0)
					System.out.println(f.format((nX/((double)nRows)*100)) + " % of file read");
				eles = line.split("\t",nCols+1);
				if(hasRowNames)
				{
					//System.out.println(eles.length+ "rowNames.length " + rowNames.length + " x=" + x);
					rowNames[x] = eles[0];
				}
				else
				{
					rowNames[x] = "Row"+Integer.toString(x);
				}
				
				for(int y =0; y < nCols-num;y++)
				{
					try
					{
						String value = "";
						if(hasRowNames && eles[y+1].length()>0)
							value= eles[y+1];
						else if (eles[y].length()>0)
							value= eles[y];
						
						values[x][y] = value;
					}catch (Exception e1)
					{
						System.out.println("FileName =  " + fileName);
						System.out.println("x =  " + x + " y = " + y + " cols[y] = " + colNames[y]);
						System.out.println("rows.length =  " + values.length + "cols.length  = " + values[x].length);
						System.out.println("eles.length =  " + eles.length);
						System.out.println("exception =  " + e1);
						values[x][y] = "";
						e1.printStackTrace();
						System.exit(666);
					}
					if(maxY > 0 && y >= maxY-1)
						break;
				}
				
				x++;
				if(maxX > 0 && x >= maxX-1)
					break;
				nX++;
			}while((line = reader.readLine())!= null);
			reader.close();
		} catch (Exception e) 
		{
			e.printStackTrace();
			if(!new File(fileName).exists())
			{
				System.out.println("File does not exist: " + fileName + "\n Exiting");
				System.exit(1);
			}
			else
				System.out.println("FileName =" + fileName);
		}

		if(verbose)
			System.out.println("Finished reading file");
		
	}
	private BufferedReader getReader(String fileName) throws IOException 
	{
		BufferedReader reader = null;
		if (fileName.endsWith(".gz"))
		{
			GZIPInputStream gzipInputStream = new GZIPInputStream(new FileInputStream(fileName));
			reader = new BufferedReader(new InputStreamReader(gzipInputStream, "US-ASCII"));
		} else 
			reader = new BufferedReader(new InputStreamReader(new FileInputStream(fileName), ENCODING), 8096);
		return reader;
	}

	public Hashtable<String,Integer> rowNamesToHash()
	{
		return namesToHash(this.rowNames);
	}
	public Hashtable<String,Integer> colNamesToHash()
	{
		return namesToHash(this.colNames);
	}
	public Hashtable<String,Integer> namesToHash(String[] names)
	{
		Hashtable<String,Integer> index = new Hashtable<String,Integer>();
		//System.out.println("names[0] = " + names[0] + " this.rowNames[0]=" + this.rowNames[0]);
		for(int x = 0; x < names.length;x++)
		{
			index.put(names[x], x);
		}
		return index;
	}
	
	
	
	public void keepRows(MatrixString m2)//Changes both matrixes if necessary, Makes sure same IDs are on same line in both files
	{
		//find all the IDs shared by both files
		Hashtable<String,Integer> allIDs = m2.rowNamesToHash();
		Hashtable<String,Integer> toKeep = new Hashtable<String,Integer>();
		int newPos = 0;
		for(int r = 0; r < this.rowNames.length;r++)
		{
			if(allIDs.containsKey(this.rowNames[r]))
			{
				toKeep.put(this.rowNames[r], newPos);
				newPos++;
			}
		}
		this.keepIDsInHash(toKeep);
		m2.keepIDsInHash(toKeep);
	}
	
	private void keepIDsInHash(Hashtable<String, Integer> toKeep) 
	{
		int total = toKeep.size();
		
		String[][] vals = new String[total][];
		String[] rN = new String[total];
		
		for(int x = 0; x < rowNames.length;x++)
		{
			if(toKeep.containsKey(rowNames[x]))
			{
				rN[toKeep.get(rowNames[x])] = rowNames[x]; 
				vals[toKeep.get(rowNames[x])] = values[x];
			}
		}
		this.values = vals;
		matrix = new GetVal(values);
		this.rowNames = rN;
	}

	public void write(String fileName)
	{
		print(rowNames.length, colNames.length,fileName,-1);
	}
	public void write(Path fileName)
	{
		print(rowNames.length, colNames.length,fileName.toString(),-1);
	}
	public void write(String fileName,int x, int y)
	{
		print(x, y, fileName,-1);
	}
	public void write(String fileName, int decimals)
	{
		print(rowNames.length, colNames.length, fileName,decimals);
	}
	public void write(Path fileName, int decimals)
	{
		print(rowNames.length, colNames.length,fileName.toString(),decimals);
	}
	public void print()
	{
		print(rowNames.length, colNames.length, null,-1);
	}
	public void print(int x, int y)
	{
		print(x, y, null,-1);
	}
	public void print(int maxX, int maxY, String fileName, int decimals)
	{
		//Timer timer = null;
	
		if(maxX > this.rowNames.length || maxX < 0)
			maxX = this.rowNames.length;
		if(maxY > this.colNames.length || maxY < 0)
			maxY = this.colNames.length;
	
		if(rowNames.length*colNames.length>10000000)//if the files is big report how far its along occasionally
		{
			//timer = new Timer();
		}
				
		try 
		{
			BufferedWriter writer = null;
			if(fileName != null)
				writer = getWriter(fileName);
			
			String format = "#";
			
			int nDecimals = decimals;
			if(nDecimals < 0)
				nDecimals = 2;
			for(int i = 0; i < nDecimals; i++)
			{
				if(i == 0) format += ".";
				format += "#";
			}
			
			DecimalFormat df = new DecimalFormat(format); 
			//DecimalFormat df = new DecimalFormat();//This causes the format to be with a , on 1,000.00, super ennoying...
			//df.setMaximumFractionDigits(2);
			for(int x = 0; x < maxX;x++)
			{
				if(x==0)
				{
					DateFormat dateFormat = new SimpleDateFormat("dd/MM/yyyy");
					Date date = new Date();
					printOrWrite(dateFormat.format(date),writer);
					
					for(int y = 0; y < maxY; y++)
					{
						printOrWrite("\t" + colNames[y],writer);
					}
					
					printOrWrite("\n",writer);
					//line = printColNames(maxY)+"\n";
				}
				printOrWrite(rowNames[x],writer);
				for(int y = 0; y < maxY; y++)
				{
					if(decimals >= 0)
					{
						printOrWrite("\t"+df.format(values[x][y]),writer);
						//line += "\t"+df.format(values[x][y]);
					}
					else
					{
						printOrWrite("\t"+values[x][y],writer);
					}
				}
				printOrWrite("\n", writer);
//				if(x%(100000000/maxY)==0 && x>0)//print time 
//				{
//					
//					double percentage = (((double)x)/((double)maxX));
//					System.out.println(df.format(percentage*100) + " % of the file has been saved");
//					timer.print(percentage);
//				}
			}
			if(writer != null)
				writer.close();
		} catch (IOException e) {System.out.println("Oh noes! An error in your print function.");e.printStackTrace();}
	}
	
	private BufferedWriter getWriter(String fileName) throws IOException {
		BufferedWriter writer = null;
		boolean gzipped = false;
		if(fileName.endsWith(".gz"))
			gzipped=true;
		
		if(fileName != null)
		{
			if (gzipped) {
				GZIPOutputStream gzipOutputStream = new GZIPOutputStream(new FileOutputStream(fileName));
				writer = new BufferedWriter(new OutputStreamWriter(gzipOutputStream), DEFAULT_BUFFER_SIZE);
			} else {
				writer = new BufferedWriter(new FileWriter(fileName), DEFAULT_BUFFER_SIZE);
			}
		}
		return writer;
	}

	private void printOrWrite(String writeString, BufferedWriter writer) throws IOException 
	{
		if(writer == null)
			System.out.print(writeString);
		else
			writer.write(writeString);
	}

	private String printColNames(int maxY) 
	{
		String line = "";
		for(int y = 0; y < maxY; y++)
		{
			line += "\t" + colNames[y];
		}
		return line;
	}
	private int getRowNumber(String fileName) 
	{
		int nLines = 0;
		try {
			BufferedReader fileReader = getReader(fileName);
			LineNumberReader lnr = new LineNumberReader(fileReader);
			lnr.skip(Long.MAX_VALUE);
			nLines = lnr.getLineNumber() + 1; //Add 1 because line index starts at 0
			fileReader.close();
		} catch (FileNotFoundException e1) {e1.printStackTrace();} catch (IOException e) {e.printStackTrace();}
		//System.out.println("lines = " + nLines);		
		return nLines;
	}

	private void revert(double[] col) 
	{
		// reverse the array
		for(int i=0;i<col.length/2;i++) {
		     // swap the elements
		     double temp = col[i];
		     col[i] = col[col.length-(i+1)];
		     col[col.length-(i+1)] = temp;
		}
		
	}
	public void transpose()
	{
		String[][] temp = new String[colNames.length][rowNames.length];
		String[] tempNames = colNames;
		colNames = rowNames;
		rowNames = tempNames;
				
		for(int x = 0; x < rowNames.length; x++)
		{
			for(int y =0; y < colNames.length; y++)
			{
				temp[x][y] = values[y][x];
			}
		}
		values = temp;
		matrix = new GetVal(values);
	}

	public void writeTransposed(String fileName) 
	{
		BufferedWriter writer = null;
		
		//if(maxX > this.rowNames.length || maxX < 0)
		int	maxX = this.rowNames.length;
		//if(maxY > this.colNames.length || maxY < 0)
		int	maxY = this.colNames.length;
		try 
		{
			if(fileName != null)
			{
				writer =  getWriter(fileName);	
			}
			String format = "#";
			
			int nDecimals = -1;
			int decimals = nDecimals;
			if(nDecimals < 0)
				nDecimals = 2;
			for(int i = 0; i < nDecimals; i++)
			{
				if(i == 0) format += ".";
				format += "#";
			}
			
			DecimalFormat df = new DecimalFormat(format); 
			//DecimalFormat df = new DecimalFormat();//This causes the format to be with a , on 1,000.00, super ennoying...
			//df.setMaximumFractionDigits(2);
			for(int x = 0; x < maxY;x++)
			{
				if(x==0)
				{
					DateFormat dateFormat = new SimpleDateFormat("dd/MM/yyyy");
					Date date = new Date();
					printOrWrite(dateFormat.format(date),writer);
					for(int y = 0; y < maxX; y++)
					{
						printOrWrite("\t" + rowNames[y],writer);
					}
					
					printOrWrite("\n",writer);
					//line = printColNames(maxY)+"\n";
				}
				printOrWrite(colNames[x],writer);
				for(int y = 0; y < maxX; y++)
				{
					if(decimals >= 0)
					{
						printOrWrite("\t"+df.format(values[y][x]),writer);
						//line += "\t"+df.format(values[x][y]);
					}
					else
					{
						printOrWrite("\t"+values[y][x],writer);
					}	
				}
				printOrWrite("\n",writer);
			}
			writer.close();
		}catch (IOException e)
		{
			e.printStackTrace();
		}	
	}

	public MatrixString getCol(int i) 
	{
		MatrixString temp = new MatrixString(this.rowNames.length,1);
		temp.rowNames = this.rowNames;
		temp.colNames = new String[]{this.colNames[i]};
		for(int r = 0; r < this.rowNames.length; r++)
		{
			temp.values[r] = this.values[r];
		}
		// TODO Auto-generated method stub
		return temp;
	}

	public MatrixString getRow(int i) 
	{
		MatrixString temp = new MatrixString(1,this.colNames.length);
		temp.rowNames = new String[]{this.rowNames[i]};
		temp.colNames = this.colNames;
		temp.values[0] = this.values[i];
		return temp;
	}
	
	public class Timer 
	{
		long startTime ;
		
		public Timer()
		{
			startTime = System.nanoTime();    
			// ... the code being measured ...    
		}
		public long runTime()
		{
			return (System.nanoTime() - startTime);
		}
		public void print()
		{
			print(0);
		}
		public void print(double percentage)
		{
			
			int seconds = (int)(runTime()/1000/1000/1000);
			int minutes = seconds/60;
			int hours = minutes/60;
			
			System.out.print("Current runtime = ");
			if(hours > 0)
				System.out.print(hours + " hours and ");
			if(minutes>0)
				System.out.print(minutes%60 + " minutes and ");
			System.out.println(seconds%60 + " seconds");
			
			if(percentage>0)
			{
				int totalRuntime = (int)(((double)seconds)/percentage);
				int timeLeft = totalRuntime-seconds;
				System.out.println("Estimated time left:" + timeLeft + " seconds. Total estimated time for this operation: " + totalRuntime);
			}	
				
			//System.out.println("Current runtime = " + runTime()/1000/1000/1000/60/24 + "days");
		}
	}

//	public void rLog(double rLog, String writeFolder, String fileName) 
//	{
//		this.logTransform(10);
//		Matrix geoMean = new Matrix(this.rowNames.length,1);
//		for(int r =0; r < this.rowNames.length; r++)
//		{
//			geoMean.values[r][0] = this.sumRow(r)/this.colNames.length;
//			
//			geoMean.values[r][0] = Math.pow(10,geoMean.values[r][0]);
//		}
//		System.out.println("GEO "+geoMean.rowNames.length);
//		if(writeFolder != null)
//			geoMean.write(writeFolder+ "geoMean.txt");
//		this.pow(10,1);//i could save the initial matrix, but this uses less memory
//					 //may be better to read in the matrix again to avoid rounding errors
//		this.divideByCol(geoMean,0);
//		
//		//need to read the matrix again here...
//		this.readFile(fileName, true, true);
//		Matrix denominators = new Matrix(this.colNames.length,1);
//
//		for(int c = 0; c < this.colNames.length; c++)
//		{
//			double[] column = new double[this.rowNames.length];
//			for(int r = 0; r < column.length; r++)
//			{
//				column[r] = (this.values[r][c]+1) / geoMean.values[r][0];// Adding +1 here just as in the log calculation earlier (if not the median can be 0 causing a division by 0 lateron)
//			}
//			Arrays.sort(column);
//			
//			denominators.values[c][0] =column[column.length/2];
//			
//			for(int r = 0; r < column.length; r++)
//			{
//				this.values[r][c] /= denominators.values[c][0];
//			}
//			
//		}
//		if(writeFolder != null)
//			denominators.write(writeFolder + "Denominators.txt");
//	}
	
	public String[] getRowHeaders()
	{
		return this.rowNames;
	}
	public Hashtable<String,Integer> getRowHash()
	{
		if(this.rowHash == null || this.rowHash.size() != this.getRowHeaders().length)
			this.rowHash = makeHash(this.rowNames);
		return this.rowHash;
	}
	public Hashtable<String,Integer> getColHash()
	{
		if(this.colHash == null || this.colHash.size() != this.getColHeaders().length)
			this.colHash = makeHash(this.colNames);
		return this.colHash;
	}
	private Hashtable<String, Integer> makeHash(String[] names) {
		Hashtable<String,Integer> hash = new Hashtable<String,Integer>();
		for(int n = 0; n < names.length; n++)
		{
			hash.put(names[n], n);
		}
		return hash;
	}
	public String[] getColHeaders()
	{
		return this.colNames;
	}
	public int rows()
	{
		return this.rowNames.length;
	}
	public int cols()
	{
		return this.colNames.length;
	}

	public void putGenesOnRows() 
	{
		if((!this.rowNames[0].contains("ENSG0") || !this.rowNames[0].contains("ENST0")) && (this.colNames[0].contains("ENSG0") || (this.colNames[0].contains("ENST0"))))
		{
			System.out.println("Putting genes on rows(transposing)");
			this.transpose();
		}
		
	}
	public void setRowHeaders(String[] rowHeaders) {
		this.rowNames = rowHeaders;
	}
	public void setColHeader(int outCol, String string) {
		this.colNames[outCol] = string;
	}
	public void set0() 
	{
		for(String[] row: values)
			for(int c = 0; c < row.length; c++)
				row[c]="0";
	}
	public String getAsLine(String[] row, String rowName) {
		StringBuilder stringLineBuilder = new StringBuilder();
		stringLineBuilder.append(rowName);
		for(String ele : row)
		{
			stringLineBuilder.append("\t");
			stringLineBuilder.append(ele);
		}
		stringLineBuilder.append("\n");
		return stringLineBuilder.toString();
	}
}
