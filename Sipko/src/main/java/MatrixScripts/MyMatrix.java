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
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.stat.descriptive.moment.Variance;

import JuhaPCA.FileUtil;
import JuhaPCA.PCA;
import no.uib.cipr.matrix.DenseMatrix;
import umcg.genetica.containers.Pair;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.Format;

public class MyMatrix
{
	//a matrix class that is not limited to max java array size

	public String firstField;
	public String[] colNames;
	public String[] rowNames;
	public double[][] values;
	public boolean verbose = false;
	private Hashtable<String, Integer> rowHash = null;
	private Hashtable<String, Integer> colHash = null;
	protected static final String ENCODING = "ISO-8859-1";
	public static final int DEFAULT_BUFFER_SIZE = 4096;
	static final Format dateFormat = new SimpleDateFormat("yyyy-MM-dd");
    static final Format timeFormat = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

	public class GetVal
	{
		double[][] values;

		GetVal(double[][] matrix)
		{
			values = matrix;
		}

		public double get(	int r,
							int c)
		{
			return values[r][c];
		}

		public void set(int r,
						int c,
						double d)
		{
			values[r][c] = d;
		};

		public void add(int r,
						int c,
						double addValue)
		{
			values[r][c] += addValue;
		}
	}

	public GetVal matrix = new GetVal(values);

	public MyMatrix()
	{

	}

	public MyMatrix(int x,
					int y)
	{
		rowNames = new String[x];
		colNames = new String[y];
		values = new double[x][y];
		matrix = new GetVal(values);

		numberNames(rowNames,
					"Row");
		numberNames(colNames,
					"Col");
	}

	public MyMatrix(int x,
					int y,
					String[] colNames,
					String[] rowNames)
	{
		this.rowNames = rowNames;
		this.colNames = colNames;
		values = new double[x][y];
		matrix = new GetVal(values);
	}

	public MyMatrix(String fileName)
	{
		readFile(	fileName,
					true,
					true);
	}

	public MyMatrix(Path fileName)
	{
		readFile(	fileName.toString(),
					true,
					true);
	}

	public MyMatrix(String fileName,
					int x,
					int y)//read file with a maximum of xRows and yCols
	{
		readFile(	fileName,
					true,
					true,
					x,
					y,0);
	}

	public MyMatrix(String fileName,
					boolean hasRowColNames)
	{
		readFile(	fileName,
					false,
					false);
	}
	public MyMatrix(String fileName,
					boolean hasRowColNames, boolean hasColNames)
	{
		readFile(	fileName,
		         	hasRowColNames,
		         	hasColNames);
	}

	public MyMatrix(String[] rowNames,
					String[] colNames,
					double[][] matrix)
	{
		this.rowNames = rowNames;
		this.colNames = colNames;
		this.values = matrix;
		this.matrix = new GetVal(values);
	}

	public MyMatrix(MatrixStruct expressionMatrixStruct)
	{
		this.rowNames = expressionMatrixStruct.getRowHeaders();
		this.colNames = expressionMatrixStruct.getColHeaders();
		this.values = new double[this.rowNames.length][];
		matrix = new GetVal(values);
		for (int r = 0; r < this.rowNames.length; r++)
			this.values[r] = expressionMatrixStruct.getRowValues(r);
	}

	double[][] DenseMatrixToValues(DenseMatrix denseMatrix)
	{
		double[][] values = new double[denseMatrix.numRows()][denseMatrix.numColumns()];
		for (int r = 0; r < denseMatrix.numRows(); r++)
			for (int c = 0; c < denseMatrix.numColumns(); c++)
				this.values[r][c] = denseMatrix.get(r,
													c);
		return values;
	}

	public MyMatrix(DenseMatrix matrix,
					String[] rowHeaders,
					String[] colHeaders)
	{
		this.setRowHeaders(rowHeaders);
		this.setColHeaders(colHeaders);
		this.values = DenseMatrixToValues(matrix);
	}

	public MyMatrix(String fileName,
					int columnToUseAsRowNames)//-1 is last column
	{
		readFile(	fileName,
					true,
					true, -1,-1,-1);
	}

	public MyMatrix(double[][] values)
	{
		this(new String [values.length],
					new String[values[0].length],
					values);
		
		this.firstField="";
		numberNames(rowNames,
				"Row");
		numberNames(colNames,
				"Col");
	}

	public void readFile(String fileName)
	{
		readFile(	fileName,
					true,
					true);
	}

	public void readFile(	String fileName,
							boolean hasRowNames,
							boolean hasColName)
	{
		readFile(	fileName,
					hasRowNames,
					hasColName,
					-1,
					-1,0);
	}

	private void numberNames(	String[] names,
								String extra)
	{
		for (int x = 0; x < names.length; x++)
		{
			names[x] = extra + "" + Integer.toString(x + 1);
		}
	}

	public void log2Transform()//add +1 before transforming
	{
		logTransform(2);
	}

	public void logTransform(int val)//add +1 before transforming
	{
		logTransform(	val,
						1);
	}

	public void logTransform(	int val,
								double addval)//add +1 before transforming
	{
		double logVal = Math.log(val);
		for (int x = 0; x < this.rowNames.length; x++)
		{
			for (int y = 0; y < this.colNames.length; y++)
			{
				this.values[x][y] = Math.log(this.values[x][y] + addval) / logVal;
			}
		}
	}

	//Matrix looks like this:
	//X n n n
	//n v v v
	//n v v v
	//n v v v
	//Where "n" is row or column name and "v" is value
	public void readFile(	String fileName,
							boolean hasRowNames,
							boolean hasColNames,
							int maxX,
							int maxY, int rowNameColumn)
	{		
		int nRows = 0;
		if (maxX < 1)
			nRows = getRowNumber(fileName);
		else
			nRows = maxX;
		
		try
		{
			
			BufferedReader reader = getReader(fileName);

			String line = null;
			line = reader.readLine();
			String[] eles = line.split("\t");
			int nCols = eles.length;
			if (nCols > maxY && maxY > 0)
				nCols = maxY;

			if (verbose)
			{
				System.out.println("Input file has the following format (rows,columns): (" + nRows + "," + nCols + ")");
				//System.out.println("This includes 1 row and 1 column for the row/col names");
			}

			String secondLine = reader.readLine();
			int secondLineCols = secondLine.split("\t").length;
			int minusCol = 1;
			if (secondLineCols == nCols + 1)//first line is missing a cell
				minusCol = 0;

			if (hasColNames && hasRowNames)
				values = new double[nRows - 1][nCols - minusCol];
			else if (hasColNames && !hasRowNames)
				values = new double[nRows-1][nCols] ;
			else if (!hasColNames && hasRowNames)
				values = new double[nRows][nCols - minusCol];
			else
				values = new double[nRows][nCols];
			
			if(rowNameColumn==-1)
				rowNameColumn=eles.length-1;
			
			matrix = new GetVal(values);
			
			colNames = new String[nCols - minusCol];
			//Deal with first line missing first cell
			int outCol=0;
			for (int y = 0; y < nCols; y++)
			{
				if(y==rowNameColumn)
				{
					if(minusCol==0)
						continue;
					firstField=eles[y];
					continue;
				}
				
				if(!hasColNames)
				{
					colNames[outCol] = "Col" + Integer.toString(y);
					firstField = null;
				}
				else
					colNames[outCol] = eles[y];
				
				outCol++;
			}
			
			if(!hasColNames) //reset reader to start of file, probably not the best way of doing it...
			{
				reader.close();
				reader = getReader(fileName);
			}
			
			if (hasColNames)
				rowNames = new String[nRows - 1];
			else
				rowNames = new String[nRows];
			
			int x = 0;
			if (hasColNames == true)
				line = secondLine;
			else
				line = reader.readLine();

			//System.out.println("");
//			double nX = 0;
//			double blockSize = 1000000000;//print after 100.000.000 cells have been read
//			int printPoint = (int) (blockSize / nCols);
//			
//			DecimalFormat f = new DecimalFormat("#.##");
			
			do
			{
				//System.out.println(line);
				if (x >= rowNames.length)// in case the last line of the file is corrupt (e.g. an incomplete tranfer) the last line will not be included
				{
					break;
				}

//				if (nX % printPoint == 0 && nX != 0)
//					System.out.println(f.format((nX / ((double) nRows) * 100)) + " % of file read");
				eles = line.split(	"\t",
									nCols + 1);
				if (hasRowNames)
				{
					//System.out.println(eles.length+ "rowNames.length " + rowNames.length + " x=" + x);
					rowNames[x] = eles[rowNameColumn];
				}
				else
				{
					rowNames[x] = "Row" + Integer.toString(x);
				}

				int yOut = 0;
				for (int y = 0; y < nCols; y++)
				{
					if(y==rowNameColumn)
						continue;
					try
					{
						double value = 0;
						if (eles[y].length() > 0)
							value = Double.parseDouble(eles[y]);
						
						values[x][yOut] = value;
					} catch (Exception e1)
					{
						System.out.println("FileName =  " + fileName);
						System.out.println("x =  " + x + " y = " + y + " cols[y] = " + colNames[y]);
						System.out.println("rows.length =  " + values.length + "cols.length  = " + values[x].length);
						System.out.println("eles.length =  " + eles.length);
						System.out.println("exception =  " + e1);
						values[x][y] = 0;
						e1.printStackTrace();
						System.exit(666);
					}
					yOut++;
					if (maxY > 0 && y >= maxY - 1)
						break;
				}

				x++;
				if (maxX > 0 && x >= maxX - 1)
					break;
//				nX++;
			} while ((line = reader.readLine()) != null);
			reader.close();
			
		} catch (Exception e)
		{
			e.printStackTrace();
			if (!new File(fileName).exists())
			{
				System.out.println("File does not exist: " + fileName + "\n Exiting");
				System.exit(1);
			}
			else
				System.out.println("FileName =" + fileName);
		}

		if (verbose)
			System.out.println("Finished reading file");
	}

	private BufferedReader getReader(String fileName) throws IOException
	{
		BufferedReader reader = null;
		if (fileName.endsWith(".gz"))
		{
			GZIPInputStream gzipInputStream = new GZIPInputStream(new FileInputStream(fileName));
			reader = new BufferedReader(new InputStreamReader(	gzipInputStream,
																"US-ASCII"));
		}
		else
			reader = new BufferedReader(new InputStreamReader(	new FileInputStream(fileName),
																ENCODING),
										8096);
		return reader;
	}

	public Hashtable<String, Integer> rowNamesToHash()
	{
		return namesToHash(this.rowNames);
	}

	public Hashtable<String, Integer> colNamesToHash()
	{
		return namesToHash(this.colNames);
	}

	public Hashtable<String, Integer> namesToHash(String[] names)
	{
		Hashtable<String, Integer> index = new Hashtable<String, Integer>();
		//System.out.println("names[0] = " + names[0] + " this.rowNames[0]=" + this.rowNames[0]);
		for (int x = 0; x < names.length; x++)
		{
			index.put(	names[x],
						x);
		}
		return index;
	}

	public void keepRows(MyMatrix m2)//Changes both matrixes if necessary, Makes sure same IDs are on same line in both files
	{
		//find all the IDs shared by both files
		Hashtable<String, Integer> allIDs = m2.rowNamesToHash();
		if(this.rowNames.length!=this.getRowHash().size() || m2.rowNames.length != allIDs.size())
		{
			System.out.println("Warning, there are multiple rows with the same IDs in either one of the matrixes. Cannot validate whether inputmatrix contains same genes/transcripts as the matrix the eigenvectors where build from");
			return;
		}
		
		Hashtable<String, Integer> toKeep = new Hashtable<String, Integer>();
		int newPos = 0;
		for (int r = 0; r < this.rowNames.length; r++)
		{
			if (allIDs.containsKey(this.rowNames[r]))
			{
				toKeep.put(	this.rowNames[r],
							newPos);
				newPos++;
			}
		}
		
		this.keepIDsInHash(toKeep);
		m2.keepIDsInHash(toKeep);
	}

	private void keepIDsInHash(Hashtable<String, Integer> toKeep)
	{
		
			int total = toKeep.size();
	
			double[][] vals = new double[total][];
			String[] rN = new String[total];
	
			for (int x = 0; x < rowNames.length; x++)
			{
				try
				{
				if (toKeep.containsKey(rowNames[x]))
				{
					rN[toKeep.get(rowNames[x])] = rowNames[x];
					vals[toKeep.get(rowNames[x])] = values[x];
				}
				}catch(Exception e)
				{
					System.out.println("x=" + x );
					System.out.println("toKeep.size()=" + toKeep.size());
					System.out.println("rowNames.length=" + rowNames.length);
					System.out.println("rN.length=" + rN.length);
					System.out.println("rowNames[x]=" + rowNames[x]);
					e.printStackTrace();
				}
			}
			this.values = vals;
			matrix = new GetVal(values);
			this.rowNames = rN;
	}

	public void write(String fileName)
	{
		print(	rowNames.length,
				colNames.length,
				fileName,
				-1);
	}

	public void write(Path fileName)
	{
		print(	rowNames.length,
				colNames.length,
				fileName.toString(),
				-1);
	}

	public void write(	String fileName,
						int x,
						int y)
	{
		print(	x,
				y,
				fileName,
				-1);
	}

	public void write(	String fileName,
						int decimals)
	{
		print(	rowNames.length,
				colNames.length,
				fileName,
				decimals);
	}

	public void write(	String fileName,
						boolean append)
	{
		print(	rowNames.length,
				colNames.length,
				fileName,
				-1,
				append);
	}
	public void write(	String fileName,
						boolean append,
						boolean writeParallel)
	{
		if(writeParallel)
		{
			ExecutorService executor = Executors.newFixedThreadPool(1);
			
			Runnable worker = new MyMatrixWorker(this, fileName);
			executor.execute(worker);
			executor.shutdown();
			System.out.println("Writing in a separate thread; continuing protocol");
		}
		else
		{
			write(fileName,append);
		}
	}

	public void write(	Path fileName,
						int decimals)
	{
		print(	rowNames.length,
				colNames.length,
				fileName.toString(),
				decimals);
	}

	public void print()
	{
		print(	rowNames.length,
				colNames.length,
				null,
				-1);
	}

	public void print(	int x,
						int y)
	{
		print(	x,
				y,
				null,
				-1);
	}

	public void print(	int maxX,
						int maxY,
						String fileName,
						int decimals)
	{
		print(	maxX,
				maxY,
				fileName,
				decimals,
				false);
	}

	public void print(	int maxX,
						int maxY,
						String fileName,
						int decimals,
						boolean append)
	{
		//Timer timer = null;

		if (maxX > this.rowNames.length || maxX < 0)
			maxX = this.rowNames.length;
		if (maxY > this.colNames.length || maxY < 0)
			maxY = this.colNames.length;

		if (rowNames.length * colNames.length > 10000000)//if the files is big report how far its along occasionally
		{
			//timer = new Timer();
		}

		try
		{
			BufferedWriter writer = null;
			if (fileName != null)
				writer = getWriter(	fileName,
									append);

			DecimalFormat df = getFormat(decimals);
			//DecimalFormat df = new DecimalFormat();//This causes the format to be with a , on 1,000.00, super ennoying...
			//df.setMaximumFractionDigits(2);
			for (int x = 0; x < maxX; x++)
			{
				if (x == 0 && !append)// write header
				{
					DateFormat dateFormat = new SimpleDateFormat("dd/MM/yyyy");
					Date date = new Date();

					if(this.firstField==null)
						this.firstField=dateFormat.format(date);
					
					printOrWrite(	this.firstField,
									writer);

					for (int y = 0; y < maxY; y++)
					{
						printOrWrite("\t",writer);
						printOrWrite(	colNames[y],
										writer);
					}

					printOrWrite(	"\n",
									writer);
					//line = printColNames(maxY)+"\n";
				}
				printOrWrite(	rowNames[x],
								writer);
				for (int y = 0; y < maxY; y++)
				{
					if (decimals >= 0)
					{
						printOrWrite("\t",writer);
						printOrWrite(	df.format(values[x][y]),
										writer);
						//line += "\t"+df.format(values[x][y]);
					}
					else
					{
						printOrWrite("\t",writer);
						printOrWrite(	Double.toString(values[x][y]),
										writer);
					}
				}
				printOrWrite(	"\n",
								writer);
				//				if(x%(100000000/maxY)==0 && x>0)//print time 
				//				{
				//					
				//					double percentage = (((double)x)/((double)maxX));
				//					System.out.println(df.format(percentage*100) + " % of the file has been saved");
				//					timer.print(percentage);
				//				}
			}
			if (writer != null)
				writer.close();
		} catch (IOException e)
		{
			System.out.println("Oh noes! An error in your print function.");
			e.printStackTrace();
		}
	}

	private DecimalFormat getFormat(int decimals)
	{
		String format = "#";

		int nDecimals = decimals;
		if (nDecimals < 0)
			nDecimals = 2;
		for (int i = 0; i < nDecimals; i++)
		{
			if (i == 0)
				format += ".";
			format += "#";
		}

		DecimalFormat df = new DecimalFormat(format);
		return df;
	}

	private BufferedWriter getWriter(String fileName) throws IOException
	{
		return getWriter(	fileName,
							false);
	}

	private BufferedWriter getWriter(	String fileName,
										boolean append) throws IOException
	{
		BufferedWriter writer = null;
		boolean gzipped = false;
		if (fileName.endsWith(".gz"))
			gzipped = true;

		if (fileName != null)
		{
			if (gzipped)
			{
				GZIPOutputStream gzipOutputStream = null;
				if (append != true)
				{
					gzipOutputStream = new GZIPOutputStream(new FileOutputStream(fileName));
				}
				else
				{
					gzipOutputStream = new GZIPOutputStream(new FileOutputStream(	fileName,
																					true));
				}
				writer = new BufferedWriter(new OutputStreamWriter(gzipOutputStream),
											DEFAULT_BUFFER_SIZE);
			}
			else
			{
				if (append != true)
				{
					writer = new BufferedWriter(new FileWriter(fileName),
												DEFAULT_BUFFER_SIZE);
				}
				else
				{
					writer = new BufferedWriter(new FileWriter(	fileName,
																true),
												DEFAULT_BUFFER_SIZE);
				}
			}
		}
		return writer;
	}

	private void printOrWrite(	String writeString,
								BufferedWriter writer) throws IOException
	{
		if (writer == null)
			System.out.print(writeString);
		else
			writer.write(writeString);
	}

	private String printColNames(int maxY)
	{
		String line = "";
		for (int y = 0; y < maxY; y++)
		{
			line += "\t" + colNames[y];
		}
		return line;
	}

	private int getRowNumber(String fileName)
	{
		int nLines = 0;
		try
		{
			BufferedReader fileReader = getReader(fileName);
			LineNumberReader lnr = new LineNumberReader(fileReader);
			lnr.skip(Long.MAX_VALUE);
			nLines = lnr.getLineNumber(); //Add 1 because line index starts at 0
			fileReader.close();
		} catch (FileNotFoundException e1)
		{
			e1.printStackTrace();
		} catch (IOException e)
		{
			e.printStackTrace();
		}
		return nLines;
	}

	public int removeNoVariance()
	{
		//important that he probes/gene/transcripts are on the X-axis (rowNames)
		ArrayList<Integer> noVarRows = new ArrayList<Integer>();
		for (int rowNumber = 0; rowNumber < this.rowNames.length; rowNumber++)//Identify all rows that have no variance
		{
			boolean hasVariance = hasVariance(rowNumber);
			if (!hasVariance)
			{
				noVarRows.add(rowNumber);
			}
		}

		if (noVarRows.size() == 0)
		{
			System.out.println("There are no rows without variance");
			return 0;
		}
		System.out.println("Removing " + noVarRows.size() + " rows without variance");
		
		double[][] remainder = new double[this.rowNames.length - noVarRows.size()][];
		String[] remainderRowNames = new String[this.rowNames.length - noVarRows.size()];

		int nextRow = 0;
		int n = 0;
		int skip = noVarRows.get(n);
		n++;
		for (int x = 0; x < this.rowNames.length; x++)
		{
			if (skip == x)
			{
				if (n < noVarRows.size())
					skip = noVarRows.get(n);
				n++;
				continue;
			}
			remainder[nextRow] = this.values[x];
			remainderRowNames[nextRow] = this.rowNames[x];
			nextRow++;
		}
		this.values = remainder;
		this.rowNames = remainderRowNames;
		matrix = new GetVal(values);
		return noVarRows.size();
	}

	private boolean hasVariance(int row)
	{
		double firstValue = this.values[row][0];
		for (int y = 0; y < this.colNames.length; y++)
		{
			if (firstValue != this.values[row][y])
				return true;
		}
		return false;
	}

	public MyMatrix calcAvgCols()
	{
		MyMatrix temp = new MyMatrix(	this.colNames.length,
										1);
		temp.rowNames = this.colNames;
		temp.colNames[0] = "Average";
		for (int y = 0; y < colNames.length; y++)
		{
			temp.values[y][0] = getAverageCol(y);
		}
		return temp;
	}

	public MyMatrix calcAvgRows()
	{
		MyMatrix temp = new MyMatrix(	this.rowNames.length,
										1);
		temp.rowNames = this.rowNames;
		temp.colNames[0] = "Average";
		for (int r = 0; r < rowNames.length; r++)
		{
			temp.values[r][0] = getAverageRow(r);
		}
		return temp;
	}

	public double getAverageRow(int r)
	{
		double avg = sumRow(r) / colNames.length;
		return avg;
	}

	public double sumRow(int r)
	{
		double sum = 0;
		for (int c = 0; c < colNames.length; c++)
		{
			sum += values[r][c];
		}
		return sum;
	}

	public MyMatrix getAverageCols()
	{
		return getAverageCols(false);
	}

	public MyMatrix getAverageCols(boolean absolute)
	{
		MyMatrix averages = new MyMatrix(	this.colNames.length,
											1);
		averages.rowNames = this.colNames;
		averages.colNames = new String[] { "ColAverages" };
		for (int c = 0; c < this.colNames.length; c++)
		{
			averages.values[c][0] = getAverageCol(	c,
													absolute);
		}
		return averages;
	}

	public double getAverageCol(int col)
	{
		return getAverageCol(	col,
								false);
	}

	public double getAverageCol(int col,
								boolean absolute)
	{
		double avg = 0;
		for (int r = 0; r < rowNames.length; r++)
		{
			if (absolute)
				avg += Math.abs(values[r][col]);
			else
				avg += values[r][col];
		}
		avg /= rowNames.length;

		return avg;
	}

	public void adjustForAverageRow(MyMatrix rowAverages,
									int row)
	{
		double average = 0;
		Hashtable<String, Integer> averagesRowHash = null;
		if (rowAverages != null)
		{
			averagesRowHash = rowAverages.getRowHash();

			average = rowAverages.values[averagesRowHash.get(this.rowNames[row])][0];
		}
		else
			average = getAverageRow(row);

		for (int c = 0; c < colNames.length; c++)
		{
			values[row][c] -= average;
		}
	}

	public void adjustForAverageAllrows(MyMatrix rowAverages) //can be called with null and will use the for the averages of the columns of this matrix
	{
		for (int r = 0; r < rowNames.length; r++)
		{
			adjustForAverageRow(rowAverages,
								r);
		}
	}

	public void adjustForAverageAllCols(MyMatrix columnAverages) //can be called with null and will use the for the averages of the columns of this matrix
	{
		for (int y = 0; y < colNames.length; y++)
		{
			adjustForAverageCol(columnAverages,
								y);
		}
	}

	public MyMatrix quantileNormVector()
	{
		MyMatrix sortedCols = new MyMatrix(	this.rowNames.length,
											this.colNames.length);
		sortedCols.rowNames = rowNames;
		sortedCols.colNames = colNames;
		for (int y = 0; y < this.colNames.length; y++)
		{
			double[] col = new double[this.rowNames.length];
			for (int x = 0; x < this.rowNames.length; x++)
			{
				col[x] = this.values[x][y];
			}
			Arrays.sort(col);

			for (int x = 0; x < this.rowNames.length; x++)
			{
				sortedCols.values[x][y] = col[x];
			}
		}
		MyMatrix qNormVector = averagesPerRow(sortedCols);
		return qNormVector;
	}

	//	private Matrix getQuantValuePerRow(Matrix sortedCols) 
	//	{
	//		Matrix qNormVector = new Matrix(sortedCols.rowNames.length, 1);
	//		for(int r = 0; r < sortedCols.rowNames.length; r++)
	//		{
	//			qNormVector.values[r][0] = getMedian(sortedCols.values[r]);
	//			//System.out.println(qNormVector.values[r][0]);
	//		}
	//		qNormVector.rowNames = sortedCols.rowNames;
	//		qNormVector.colNames = new String[]{"Averages"};
	//		return qNormVector;
	//	}
	//
	//	private double getMedian(double[] row) 
	//	{
	//		double median = 0;
	//		Arrays.sort(row);
	//		median = row[row.length/2];
	//		return median;
	//	}

	public MyMatrix quantileNormVectorCenteredLog2(String saveNameQnormVector)//Calculated the averages for each rank
	{
		MyMatrix qNormVector = quantileNormVector();
		qNormVector.log2Transform();
		qNormVector.adjustForAverageAllCols(null);
		qNormVector.write(saveNameQnormVector.replace(	".txt",
														"_TransformedCentered.txt"));
		System.out.println("Quantile normalization vector saved to " + saveNameQnormVector);
		qNormVector.rowNames = this.rowNames;//just need these to determine which genes to include in quantile normalization.
		qNormVector.colNames[0] = "centered log2 transformed quantile vector. RowNames have nothing to do with the quant vector!";
		return qNormVector;
	}

	private void revert(double[] col)
	{
		// reverse the array
		for (int i = 0; i < col.length / 2; i++)
		{
			// swap the elements
			double temp = col[i];
			col[i] = col[col.length - (i + 1)];
			col[col.length - (i + 1)] = temp;
		}

	}

	public MyMatrix averagesPerRow(MyMatrix sortedCols)
	{
		return this.calcAvgRows();
	}

	public void quantileNormAdjust(MyMatrix qNormVector)//this comes after calculating the determining the means for each rank
	{
		if (qNormVector.rowNames.length != this.rowNames.length)
		{
			System.out.println("Dimensions incorrect, transposing");
			this.transpose();
			if (qNormVector.rowNames.length != this.rowNames.length)
			{
				System.out.println("Dimensions still incorrect incorrect: \n" + " this.rowNames.length = " + this.rowNames.length + " this.colNames.length " + this.rowNames.length + "/n" + " qNormVector.rowNames.length = " + qNormVector.rowNames.length);
			}
		}

		double[][] column = new double[this.rowNames.length][2];
		for (int y = 0; y < this.colNames.length; y++)//sort each column
		{
			for (int x = 0; x < this.rowNames.length; x++)
			{
				column[x][0] = this.values[x][y];
				column[x][1] = x;
				//System.out.println(this.rowNames[x] +" "+ x + " val:" +this.values[x][y]);
			}

			Arrays.sort(column,
						new Comparator<double[]>()
						{
							public int compare(	double[] s1,
												double[] s2)
							{
								if (s1[0] > s2[0])
									return 1; // tells Arrays.sort() that s1 comes after s2
								else if (s1[0] < s2[0])
									return -1; // tells Arrays.sort() that s1 comes before s2
								else
								{
									return 0;
								}
							}
						});

			int qNormVectorIndAdj = 0;//to adjust for values that are equal to another value in the same column
			int pos = 0;
			double val = 0;
			for (int x = 0; x < column.length; x++)
			{
				//far to complicated if else to achieve something simple, but its late :S	
				if (x > 0 && column[x][0] != column[x - 1][0] && (x + 1 < column.length && column[x][0] != column[x + 1][0]))//if this value is not equal to next or previous value
				{
					qNormVectorIndAdj = 0;
					pos = x;
					val = qNormVector.values[pos][0];
				}
				else
				{
					if (x == 0 || (x > 0 && column[x][0] != column[x - 1][0]))//if this value is different from the last one recalculate the number to use (or if it is the first value)
					{
						qNormVectorIndAdj = 0;
						int extra = 1;
						while (x + extra < column.length && column[x][0] == column[x + extra][0])
						{
							qNormVectorIndAdj++;
							extra++;
						}
						pos = x + (qNormVectorIndAdj / 2);
					}

					if (qNormVectorIndAdj % 2 == 0 || pos == qNormVector.values.length - 1)
						val = qNormVector.values[pos][0];
					else
					{
						val = (qNormVector.values[pos][0] + qNormVector.values[pos + 1][0]) / 2;
					}
				}

				//System.out.println("x =" + x + "Initial value: " + column[x][0] + "vNew value: " + val + " outputRow: " + (int) column[x][1] + " quantVectorRow: " + x + " c: " + y + " rowname = " + this.rowNames[(int) column[x][1]]);
				this.values[(int) column[x][1]][y] = val;
			}
		}

	}

	public void transform(MyMatrix mat)//transforms the values of this matrix to the new PC space. not to be confused with transpose ;)
	{
		//If the expression matrix has the samples on the Y-axis and the genes on the X-axis
		//If the vectors of the different PCs are in columns (where each factor is the PC for each gene)
		//multiply each column of the expression matrix with each row of the PC matrix
		MyMatrix temp = new MyMatrix(	mat.rowNames.length,
										this.colNames.length);//rowNames should be the samples
		temp.colNames = this.colNames;
		//give the rowNames PC1, PC2, PC3, etc..
		for (int x = 0; x < temp.rowNames.length; x++)
		{
			temp.rowNames[x] = "PC" + Integer.toString(x + 1);
		}
		Timer timer = new Timer();
		//Calculate the new cooridinate for 1 sample each iteration
		for (int y = 0; y < this.colNames.length; y++)
		{
			if (y % 1000 == 0)//(this.colNames.length/10000)
			{
				double percentage = ((double) (y) / ((double) (this.colNames.length)));
				System.out.println(((int) (percentage * 100)) + "% completed");
				timer.print(percentage);
			}
			for (int x = 0; x < mat.rowNames.length; x++)
			{
				//System.out.println("x = " + x + "/" + mat.rowNames.length);
				temp.values[x][y] = multiplyVectors(mat,
													x,
													y);
			}
		}

		this.values = temp.values;
		this.rowNames = temp.rowNames;

		//resulting file (temp), has PCs on RowNames, samples as colNames
	}

	private double multiplyVectors(	MyMatrix mat,
									int col,
									int y2) //multiply a column of the expression data with a row of the PC matrix (mat= PC matrix)
	{
		double result = 0;
		double[] product = new double[mat.colNames.length];
		for (int y = 0; y < mat.colNames.length; y++)
		{
			//System.out.println("col = " + col + " y = " + y+ " this.values[y][col] " + this.values[y][y2] + " y2 = " + y2);
			product[y] = this.values[y][y2] * mat.values[col][y];//x in values is the columns	
		}
		double sum = sumArray(product);
		result = sum;
		return result;
	}

	private double sumArray(double[] numbers)
	{
		double sum = 0;
		for (int x = 0; x < numbers.length; x++)
		{
			sum += numbers[x];
		}
		return sum;
	}

	public void transpose()
	{
		double[][] temp = new double[colNames.length][rowNames.length];
		String[] tempNames = colNames;
		colNames = rowNames;
		rowNames = tempNames;

		for (int x = 0; x < rowNames.length; x++)
		{
			for (int y = 0; y < colNames.length; y++)
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
		int maxX = this.rowNames.length;
		//if(maxY > this.colNames.length || maxY < 0)
		int maxY = this.colNames.length;
		try
		{
			if (fileName != null)
			{
				writer = getWriter(fileName);
			}
			String format = "#";

			int nDecimals = -1;
			int decimals = nDecimals;
			if (nDecimals < 0)
				nDecimals = 2;
			for (int i = 0; i < nDecimals; i++)
			{
				if (i == 0)
					format += ".";
				format += "#";
			}

			DecimalFormat df = new DecimalFormat(format);
			//DecimalFormat df = new DecimalFormat();//This causes the format to be with a , on 1,000.00, super ennoying...
			//df.setMaximumFractionDigits(2);
			for (int x = 0; x < maxY; x++)
			{
				if (x == 0)
				{
					DateFormat dateFormat = new SimpleDateFormat("dd/MM/yyyy");
					Date date = new Date();
					printOrWrite(	dateFormat.format(date),
									writer);
					for (int y = 0; y < maxX; y++)
					{
						printOrWrite("\t",writer);
						printOrWrite(	rowNames[y],
										writer);
					}

					printOrWrite(	"\n",
									writer);
					//line = printColNames(maxY)+"\n";
				}
				printOrWrite(	colNames[x],
								writer);
				for (int y = 0; y < maxX; y++)
				{
					if (decimals >= 0)
					{
						printOrWrite("\t",writer);
						printOrWrite(	df.format(values[y][x]),
										writer);
						//line += "\t"+df.format(values[x][y]);
					}
					else
					{
						printOrWrite("\t",writer);
						printOrWrite(	Double.toString(values[y][x]),
										writer);
					}
				}
				printOrWrite(	"\n",
								writer);
			}
			writer.close();
		} catch (IOException e)
		{
			e.printStackTrace();
		}
	}

	public MyMatrix getCol(int i)
	{
		MyMatrix temp = new MyMatrix(	this.rowNames.length,
										1);
		temp.rowNames = this.rowNames;
		temp.colNames = new String[] { this.colNames[i] };
		for (int r = 0; r < this.rowNames.length; r++)
		{
			temp.values[r] = this.values[r];
		}
		// TODO Auto-generated method stub
		return temp;
	}

	public MyMatrix getRow(int i)
	{
		MyMatrix temp = new MyMatrix(	1,
										this.colNames.length);
		temp.rowNames = new String[] { this.rowNames[i] };
		temp.colNames = this.colNames;
		temp.values[0] = this.values[i];
		return temp;
	}

	public class Timer
	{
		long startTime;

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

			int seconds = (int) (runTime() / 1000 / 1000 / 1000);
			int minutes = seconds / 60;
			int hours = minutes / 60;

			System.out.print("Current runtime = ");
			if (hours > 0)
				System.out.print(hours + " hours and ");
			if (minutes > 0)
				System.out.print(minutes % 60 + " minutes and ");
			System.out.println(seconds % 60 + " seconds");

			if (percentage > 0)
			{
				int totalRuntime = (int) (((double) seconds) / percentage);
				int timeLeft = totalRuntime - seconds;
				System.out.println("Estimated time left:" + timeLeft + " seconds. Total estimated time for this operation: " + totalRuntime);
			}

			//System.out.println("Current runtime = " + runTime()/1000/1000/1000/60/24 + "days");
		}
	}

	public void correctForTotalReadCount(double endCounts)
	{
		//according to http://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29 (voom)
		//log 2 transformation is done after/outside this function though
		for (int c = 0; c < this.colNames.length; c++)
		{
			double total = 0;
			for (int r = 0; r < this.rowNames.length; r++)
			{
				total += this.values[r][c];
			}
			total++;
			for (int r = 0; r < this.rowNames.length; r++)
			{
				this.values[r][c] = (this.values[r][c] + 0.5) / total * endCounts;//(geneExpression+0.5)/totalExpresion*10^6
			}
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

	public Hashtable<String, Integer> getRowHash()
	{
		if (this.rowHash == null || this.rowHash.size() != this.getRowHeaders().length)
			this.rowHash = makeHash(this.rowNames);
		return this.rowHash;
	}

	public Hashtable<String, Integer> getColHash()
	{
		if (this.colHash == null || this.colHash.size() != this.getColHeaders().length)
			this.colHash = makeHash(this.colNames);
		return this.colHash;
	}

	static public Hashtable<String, Integer> makeHash(String[] names)
	{
		Hashtable<String, Integer> hash = new Hashtable<String, Integer>();
		for (int n = 0; n < names.length; n++)
		{
			if (hash.containsKey(names[n]))
				System.out.println("Warning duplicate rowNames are found; Rowhash may refer to incorrect indexes\t" + n);
			if (hash.containsKey(""))
				System.out.println("Warning EMPTY rowNames are found!\t");
			hash.put(	names[n],
						n);
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

	private void divideByCol(	MyMatrix denominators,
								int col)
	{
		for (int c = 0; c < this.colNames.length; c++)
		{
			for (int r = 0; r < this.rowNames.length; r++)
			{
				this.values[r][c] /= denominators.values[c][col];
			}
		}
	}

	private void pow(int i)
	{
		pow(i,
			0);
	}

	private void pow(	int i,
						int minus)
	{
		for (int x = 0; x < this.rowNames.length; x++)
		{
			for (int y = 0; y < this.colNames.length; y++)
			{
				this.values[x][y] = Math.pow(	10,
												this.values[x][y])
						- minus;
			}
		}
	}

	public void setRowHeaders(String[] rowHeaders)
	{
		this.rowNames = rowHeaders;
	}

	public void setColHeaders(String[] colHeaders)
	{
		this.colNames = colHeaders;
	}

	public void roundValues()
	{
		for (int x = 0; x < this.rows(); x++)
		{
			for (int y = 0; y < this.cols(); y++)
				this.matrix.set(x,
								y,
								(double) (int) this.matrix.get(	x,
																y));
		}
	}

	public void expressionToRank(MyMatrix vector)
	{
		expressionToRank(	vector,
							0);
	}

	public void expressionToRank(	MyMatrix vector,
									double min)
	{
		double[][] column = new double[this.rows()][2];
		for (int c = 0; c < this.cols(); c++)//sort each column
		{
			for (int x = 0; x < this.rows(); x++)
			{
				column[x][0] = this.matrix.get(	x,
												c);
				column[x][1] = x;
			}

			Arrays.sort(column,
						new Comparator<double[]>()
						{
							public int compare(	double[] s1,
												double[] s2)
							{
								if (s1[0] > s2[0])
									return 1; // tells Arrays.sort() that s1 comes after s2
								else if (s1[0] < s2[0])
									return -1; // tells Arrays.sort() that s1 comes before s2
								else
								{
									return 0;
								}
							}
						});

			int qNormVectorIndAdj = 0;//to adjust for values that are equal to another value in the same column
			int pos = 0;
			double val = 0;
			if (vector == null)
			{
				vector = new MyMatrix(	this.rows(),
										1);
				for (int r = 0; r < vector.rows(); r++)
					vector.matrix.set(	r,
										0,
										r);
			}

			for (int r = 0; r < column.length; r++)
			{
				//far to complicated if else to achieve something simple, but its late :S	
				if (r > 0 && column[r][0] != column[r - 1][0] && (r + 1 < column.length && column[r][0] != column[r + 1][0]))//if this value is not equal to next or previous value
				{
					qNormVectorIndAdj = 0;
					pos = r;
					val = vector.matrix.get(pos,
											0);
				}
				else
				{
					if (r == 0 || (r > 0 && column[r][0] != column[r - 1][0]))//if this value is different from the last one recalculate the number to use (or if it is the first value)
					{
						qNormVectorIndAdj = 0;
						int extra = 1;
						while (r + extra < column.length && column[r][0] == column[r + extra][0])
						{
							qNormVectorIndAdj++;
							extra++;
						}
						pos = r + (qNormVectorIndAdj / 2);
					}

					if (qNormVectorIndAdj % 2 == 0 || pos == vector.rows() - 1)
						val = vector.matrix.get(pos,
												0);
					else
					{
						val = (vector.matrix.get(	pos,
													0)
								+ vector.matrix.get(pos + 1,
													0))
								/ 2;
					}
				}

				//replace the actual value
				double before = this.matrix.get((int) column[r][1],
												c);
				if (min != 0 && before < min)
					this.matrix.set((int) column[r][1],
									c,
									0);
				else
					this.matrix.set((int) column[r][1],
									c,
									val);
				double after = this.matrix.get(	(int) column[r][1],
												c);
			}
		}
	}

	private String[] makeHeaders(	String prefix,
									int nCols)
	{
		String[] headers = new String[nCols];
		for (int h = 0; h < headers.length; h++)
			headers[h] = prefix + Integer.toString(h + 1);
		return headers;
	}

	//	private void setMatrix(double[][] expression)
	//	{
	//		matrix = new DenseMatrix(	expression.length,
	//									expression[0].length);
	//		for (int r = 0; r < expression.length; r++)
	//		{
	//			for (int c = 0; c < expression[0].length; c++)
	//			{
	//				matrix.set(	r,
	//							c,
	//							expression[r][c]);
	//			}
	//		}
	//	}

	public void addAveragesCols(MatrixStruct colAverages)
	{
		for (int c = 0; c < this.getColHeaders().length; c++)
		{
			double average = colAverages.matrix.get(c,
													0);
			for (int r = 0; r < this.getRowHeaders().length; r++)
			{
				this.matrix.add(r,
								c,
								average);
			}
		}
	}

	public void addAveragesRows(MatrixStruct rowAverages)
	{
		for (int r = 0; r < this.getRowHeaders().length; r++)
		{
			double average = rowAverages.matrix.get(rowAverages.getRowHash().get(this.getRowHeaders()[r]),
													0);
			for (int c = 0; c < this.getColHeaders().length; c++)
			{
				//System.out.println("" + rowAverages.rowHeaders[rowAverages.rowHash.get(this.rowHeaders[r])] + " "+ this.rowHeaders[r] + " average = " + average);
				this.matrix.add(r,
								c,
								average);
			}
		}
	}

	public int countCol(int col,
						double cutoff)
	{
		int n = 0;
		for (int r = 0; r < this.rows(); r++)
		{
			if (this.matrix.get(r,
								col) > cutoff)
				n++;
		}
		return n;
	}

	private MatrixStruct getSubMatrix(	int rows,
										int cols)
	{
		if (rows == -1 || this.rows() < rows)
			rows = this.rows();
		if (cols == -1 || this.cols() < cols)
			cols = this.cols();
		MatrixStruct subMatrix = new MatrixStruct(	rows,
													cols);
		for (int r = 0; r < subMatrix.rows(); r++)
		{
			subMatrix.setRowHeader(	r,
									this.getRowHeaders()[r]);
			for (int c = 0; c < subMatrix.cols(); c++)
			{
				if (r == 0)
					subMatrix.setColHeader(	c,
											this.getColHeaders()[c]);
				subMatrix.matrix.set(	r,
										c,
										this.matrix.get(r,
														c));
			}
		}
		return subMatrix;
	}

	public void setZeroRows(ArrayList<Integer> PCsToCorrect)
	{
		for (int pc : PCsToCorrect)
		{
			for (int c = 0; c < this.getColHeaders().length; c++)
				this.matrix.set(pc,
								c,
								0);
		}

	}

	public void setZeroRows(String[] PCsToCorrect)
	{
		for (String PC : PCsToCorrect)
		{
			if (this.getRowHash().get(PC) == null)
			{
				System.out.println("PC not found: " + PC);
				continue;
			}
			int r = this.getRowHash().get(PC);
			System.out.println(r);
			for (int c = 0; c < this.getColHeaders().length; c++)
				this.matrix.set(r,
								c,
								0);
		}
	}

	public double[] getRowValues(int i)
	{
		double[] row = new double[getColHeaders().length];
		for (int c = 0; c < getColHeaders().length; c++)
		{
			row[c] = matrix.get(i,
								c);
		}
		return row;
	}

	public MyMatrix copy()
	{
		MyMatrix copy = new MyMatrix();
		copy.setRowHeaders(deepStringCopy(this.getRowHeaders()));
		copy.setColHeaders(deepStringCopy(this.getColHeaders()));
		//deepcopy
		copy.values = new double[this.rows()][this.cols()];
		copy.matrix = new GetVal(copy.values);
		for (int r = 0; r < this.rows(); r++)
		{
			copy.setRow(r,
						this.getRowHeaders()[r],
						this.getRowValues(r));
		}
		return copy;
	}

	private String[] deepStringCopy(String[] headers2)
	{
		String[] copy = new String[headers2.length];
		for (int e = 0; e < headers2.length; e++)
		{
			copy[e] = headers2[e];
		}
		// TODO Auto-generated method stub
		return copy;
	}

	public void setRow(	int rowNumber,
						String rowName,
						double[] rowValues)
	{
		this.getRowHeaders()[rowNumber] = rowName;
		this.rowHash = null;
		//this is bad as it remakes the whole has everytime you change a row:
		//this more then doubles the time this function takes
		
		this.setRowValues(rowNumber ,rowValues);
	}
	

	private void setRowValues(int rowNumber, double[] rowValues)
	{
		for (int c = 0; c < rowValues.length; c++)
		{
			this.matrix.set(rowNumber,
							c,
							rowValues[c]);
		}// TODO Auto-generated method stub
		
	}

	public void keepRows1Matrix(MyMatrix sample)
	{
		//keeps only the IDs in Sample
		HashMap<String, Integer> toKeep = new HashMap<String, Integer>();
		int newPos = 0;
		for (int r = 0; r < this.rows(); r++)
		{
			if (sample.getRowHash().containsKey(this.getRowHeaders()[r]))
			{
				toKeep.put(	this.getRowHeaders()[r],
							newPos);
				newPos++;
			}
		}
		//remove the redundant IDs from the matrix
		this.keepIDs(toKeep);
	}

	public void keepRows(MatrixStruct sample)
	{
		//the order that is preserved is the one from the matrix that calls the function.
		//Does not work if either of the matrixes has multiple rows with the same IDs
		HashMap<String, Integer> toKeep = new HashMap<String, Integer>();
		int newPos = 0;
		for (int r = 0; r < this.rows(); r++)
		{
			if (sample.getRowHash().containsKey(this.getRowHeaders()[r]))
			{
				toKeep.put(	this.getRowHeaders()[r],
							newPos);
				newPos++;
			}
		}
		//remove the redundant IDs from the matrix
		this.keepIDs(toKeep);
		sample.keepIDs(toKeep);
	}

	//also puts the output matrix rows in the same order
	public void keepIDs(HashMap<String, Integer> toKeep)
	{
		MyMatrix temp = new MyMatrix(	toKeep.size(),
										this.cols());
		for (int r = 0; r < this.rows(); r++)
		{
			log("rowHeaders" + this.getRowHeaders()[r] + "\t" + toKeep.get(this.getRowHeaders()[r]));
			if (toKeep.get(this.getRowHeaders()[r]) == null)
				continue;
			
			int row = toKeep.get(this.getRowHeaders()[r]);
			
			temp.setRow(row,
						this.getRowHeaders()[r],
						this.getRowValues(r));
		}
		this.setRowHeaders(temp.getRowHeaders());
		this.matrix = temp.matrix;
	}

	public MyMatrix stDevRows()
	{
		MyMatrix stDev = new MyMatrix(	this.rows(),
										1);
		stDev.setColHeaders(new String[] { "stDev" });
		stDev.setRowHeaders(this.getRowHeaders());

		for (int r = 0; r < this.rows(); r++)
		{
			double[] rowValues = this.getRowValues(r);
			//mean = org.apache.commons.math3.stat.StatUtils.mean(rowValues);
			double variance = org.apache.commons.math3.stat.StatUtils.variance(rowValues);
			double standardDev = java.lang.Math.pow(variance,
													0.5);
			stDev.matrix.set(	r,
								0,
								standardDev);
		}
		return stDev;
	}

	public MatrixStruct stDevCols()
	{
		MatrixStruct stDev = new MatrixStruct(	this.cols(),
												1);
		stDev.setColHeaders(new String[] { "stDevCols" });
		stDev.setRowHeaders(this.getColHeaders());

		for (int c = 0; c < this.cols(); c++)
		{
			double[] values = this.getColValues(c);
			//mean = org.apache.commons.math3.stat.StatUtils.mean(rowValues);
			double variance = org.apache.commons.math3.stat.StatUtils.variance(values);
			double standardDev = java.lang.Math.pow(variance,
													0.5);
			stDev.matrix.set(	c,
								0,
								standardDev);
		}
		return stDev;
	}
	public double[] getColValues(String colName)
	{
		if(!this.getColHash().containsKey(colName))
			return null;
		
		int c = this.getColHash().get(colName);
		return this.getColValues(c);
	}

	public Double[] getColValues(String colName, boolean nonPrimitives)
	{
		if(!this.getColHash().containsKey(colName))
			return null;
		
		int c = this.getColHash().get(colName);
		return this.getColValues(c, true);
	}
	public Double[] getColValues(int c, boolean nonPrimitives)
	{
		Double[] col = new Double[this.rows()];

		for (int r = 0; r < this.rows(); r++)
		{
			col[r] = this.matrix.get(	r,
										c);
		}
		return col;
	}
	
	public double[] getColValues(int c)
	{
		double[] col = new double[this.rows()];

		for (int r = 0; r < this.rows(); r++)
		{
			col[r] = this.matrix.get(	r,
										c);
		}
		return col;
	}

	public void divideBy(	MyMatrix denoms,
							boolean rows)
	{
		if (denoms.rows() != this.rows() && rows)
		{
			System.out.println("Warning! These matrixes do not have equal row lengths: rows input: " + this.rows() + "rows denominators:" + denoms.rows());
			System.out.println("Aligning rowNames");
			denoms.keepRows(this);
		}
		if (denoms.rows() != this.cols() && !rows)
		{
			System.out.println("Warning! These rows of the denominator vector is not of the same length as the cols of the matrix to adjust. column input: " + this.cols() + "rows denominators:" + denoms.rows());
		}

		//I should use vector multiplications Vector no.uib.cipr.matrix.Matrix.mult(Vector x, Vector y)
		if (denoms.cols() == 1)//divide every column by the single column in denoms
		{
			for (int r = 0; r < this.rows(); r++)
			{
				for (int c = 0; c < this.cols(); c++)
				{
					double currentVal = this.matrix.get(r,
														c);
					if (rows)
					{
						this.matrix.set(r,
										c,
										currentVal / denoms.matrix.get(	r,
																		0));
					}
					else
					{
						this.matrix.set(r,
										c,
										currentVal / denoms.matrix.get(	c,
																		0));
					}
				}
			}
		}
	}

	public void setCol(	double[] colValues,
						int c)
	{
		for (int r = 0; r < this.rows(); r++)
		{
			this.matrix.set(r,
							c,
							colValues[r]);
		}
	}

	public double[][] getMatrix()
	{
		double[][] matrix = new double[this.rows()][this.cols()];
		for (int r = 0; r < this.rows(); r++)
		{
			matrix[r] = this.getRowValues(r);
		}
		return matrix;
	}

	public void sortCol(int col)
	{
		sortCol(col,
				1);
	}

	public void sortCol(int col,
						int opposite)
	{
		class Row
		{
			String rowName;
			double[] values;

			Row(String rN,
				double[] vals)
			{
				rowName = rN;
				values = vals;
			}
		}
		ArrayList<Row> rows = new ArrayList<Row>();
		for (int r = 0; r < this.rows(); r++)
		{
			rows.add(new Row(	this.getRowHeaders()[r],
								this.getRowValues(r)));
		}

		Collections.sort(	rows,
							new Comparator<Row>()
							{
								@Override
								public int compare(	Row r1,
													Row r2)//pass this function as an argument
								{
									if (r1.values[col] > r2.values[col])
										return -1 * opposite;
									else if (r1.values[col] < r2.values[col])
										return 1 * opposite;
									else
										return 0;
								}
							});
		//rows are now sorted on cols, put it back in the actual matrix
		for (int r = 0; r < this.rows(); r++)
		{
			this.setRow(r,
						rows.get(r).rowName,
						rows.get(r).values);
		}
		this.rowHash = makeHash(this.getRowHeaders());
	}

	public MyMatrix mergeColumns(MyMatrix addition)
	{
		return mergeColumns(addition,
							false);
	}

	public MyMatrix mergeColumns(	MyMatrix addition,
									boolean keepAll)//adds columns to a matrix
	{
		if(!keepAll)
			return mergeColumns(addition, 0);
		else //keep all from this matrix, but not from addtion
			return mergeColumns(addition, 1);
	}
	
	public MyMatrix mergeColumns(	MyMatrix addition,
									int keepAll)//adds columns to a matrix
	{
		int n = 0;
		Hashtable<String,Integer> newRowHash = this.getRowHash();
		if (keepAll == 2)//only keeps all form this matrix, not from the addition matrix
		{
			newRowHash.putAll(addition.getRowHash());
			n = newRowHash.size();
		}
		else if(keepAll == 0)
		{
			for (int r = 0; r < this.rows(); r++)
			{
				if (!addition.getRowHash().containsKey(this.getRowHeaders()[r]))//if this row does not exist in the other file skip it
				{
					newRowHash.remove(this.getRowHeaders()[r]);
				}
			}
		}

		//create the new matrix
		MyMatrix result = new MyMatrix(	newRowHash.size(),
										this.cols() + addition.cols());
		
		//set the column headers
		int outC = 0;
		for (int c = 0; c < this.cols(); c++)//copy original column headers
		{
			result.setColHeader(outC,
								this.getColHeaders()[c]);
			outC++;
		}
		for (int c = 0; c < addition.cols(); c++)//add addtional column headers
		{
			result.setColHeader(outC,
								addition.getColHeaders()[c]);
			outC++;
		}
		result.colHash = MatrixStruct.makeHash(result.getColHeaders());
		//set the rowHeaders
		int outRow=0;
		Enumeration<String> rowNames = newRowHash.keys();
		
		int r = 0;
		while(rowNames.hasMoreElements())
		{
			result.rowNames[r]=rowNames.nextElement();
			r++;
		}
			
		result.rowHash=MyMatrix.makeHash(result.rowNames);
		
		//add the new rows
		for (r = 0; r < result.rows(); r++)
		{
			double[] row = new double[this.cols() + addition.cols()];
			
			row = addRowValuesToNewRow(this,row,0,result.getRowHeaders()[r]);
			row = addRowValuesToNewRow(addition,row,this.cols(),result.getRowHeaders()[r]);//add them at the end of what was added already
			
			result.setRowValues(r,
							row);
		}
		System.out.println("All colulmns merged and " + (this.rows() - result.rows()) + " rows added compared to first input file (or removed if negative)");
		return result;
	}

	private double[] addRowValuesToNewRow(	MyMatrix myMatrix,
											double[] row, int newColStart, String rowHeader)
	{
		if(myMatrix.getRowHash().containsKey(rowHeader))//if it exists in the added file
		{
			int r=myMatrix.getRowHash().get(rowHeader);
			for (int c = 0; c < myMatrix.cols(); c++)
			{
				row[c+newColStart] = myMatrix.matrix.get(r,
											c);
			}
		}
		return row;
	}

	public void setColHeader(	int c,
								String name)
	{
		this.colNames[c] = name;
		this.colHash = null;
	}

	public void setRowHeader(	int r,
								String name)
	{
		this.rowNames[r] = name;
		this.rowHash = null;
	}

	//	public void TPM(MatrixStruct geneLengths, boolean rows)//assumes geneLengths has the genes in the same order as this matrix
	//	{
	//		if(rows)//samples are on rows instead of columns
	//		{
	//			//correct for gene length
	//			this.divideBy(geneLengths, false);
	//			for(int r = 0; r < this.rows(); r++)
	//			{
	//				double total = 0;
	//				for(int c = 0; c < this.cols();c++)
	//				{
	//					total+= this.matrix.get(r, c);
	//				}
	//				if(total == 0)
	//				{
	//					System.out.println("Warning: This sample has no expression");
	//					return;
	//				}
	//				double correction = total/1000000;
	//				for(int c = 0; c < this.cols();c++)
	//				{
	//					this.matrix.set(r,c,this.matrix.get(r, c)/correction);
	//				}
	//				
	//			}
	//		}
	//		else//samples are on columns
	//		{
	//			
	//		}
	//		
	//	}

	public MyMatrix getAveragesPerRow()
	{
		return getAveragesPerRow(this);
	}

	public MyMatrix getAveragesPerRow(MyMatrix data)
	{
		MyMatrix qNormVector = new MyMatrix(data.rows(),
											1);
		for (int x = 0; x < data.rows(); x++)
		{
			double average = 0;
			for (int y = 0; y < data.cols(); y++)
			{
				average += data.matrix.get(	x,
											y);
			}
			average /= data.cols();

			qNormVector.matrix.set(	x,
									0,
									average);
		}
		qNormVector.setRowHeaders(data.getRowHeaders());
		qNormVector.setColHeaders(new String[] { "Averages" });
		return qNormVector;
	}

	public MyMatrix getAveragesPerCol()
	{
		return getAveragesPerCol(this);
	}

	public MyMatrix getAveragesPerCol(MyMatrix data)
	{
		MyMatrix temp = new MyMatrix(	data.cols(),
										1);
		temp.setRowHeaders(data.getColHeaders());
		temp.setColHeader(	0,
							"Average");
		for (int c = 0; c < data.cols(); c++)
		{
			temp.matrix.set(c,
							0,
							data.getAverageCol(c));
		}
		return temp;
	}

	public void adjustForAverageAllGenes(MyMatrix rowAverages) //can be called with null and will use the for the averages of the columns of this matrix
	{
		for (int r = 0; r < rows(); r++)
		{
			adjustForAverageRow(rowAverages,
								r);
		}
	}

	public void adjustForAverageAllSamples(MyMatrix columnAverages) //can be called with null and will use the for the averages of the columns of this MatrixStruct
	{
		for (int y = 0; y < cols(); y++)
		{
			adjustForAverageCol(columnAverages,
								y);
		}
	}

	public void adjustForAverageCol(MyMatrix columnAverages,
									int col)
	{
		double average = 0;
		if (columnAverages == null)
			average = getAverageCol(col);
		else
			average = columnAverages.matrix.get(col,
												0);
		for (int x = 0; x < rows(); x++)
		{
			//			if(col <= 1 && x < 10)
			//				System.out.println("rowname" + this.getRowHeaders()[x] + " colname= "+ this.getColHeaders()[col] + " this.matrix.get(x,col)=" + this.matrix.get(x,col) + " average =" + average);
			double value = this.matrix.get(	x,
											col);
			if (Double.isFinite(value))
				this.matrix.set(x,
								col,
								value - average);
			else
				this.matrix.set(x,
								col,
								0);
		}
	}

	public void adjustForAverageRow(MatrixStruct rowAverages,
									int row)
	{
		double average = 0;
		Hashtable<String, Integer> averagesRowHash = null;
		if (rowAverages != null)
		{
			averagesRowHash = rowAverages.getRowHash();
			average = rowAverages.matrix.get(	averagesRowHash.get(this.getRowHeaders()[row]),
												0);
		}
		else
			average = getAverageRow(row);
		for (int c = 0; c < this.cols(); c++)
		{
			double value = this.matrix.get(	row,
											c);
			if (Double.isFinite(value))
				this.matrix.set(row,
								c,
								value - average);
			else
				this.matrix.set(row,
								c,
								0);
		}
	}

	public void removeLowExpression(String removedGenesFN,
									double minSamplesExpressed,
									double minExpression) throws IOException
	{
		//important that he probes/gene/transcripts are on the rowNames
		if (minSamplesExpressed < 1)
			minSamplesExpressed = this.cols();
		ArrayList<Integer> noVarRows = new ArrayList<Integer>();//contains rows that should be removed
		for (int x = 0; x < this.rows(); x++)//Identify all rows that have no variance or expressed in less then "minSamplesExpressed" samples
		{
			int[] res = hasVariance(x,
									minExpression);
			int countExpressed = res[1];
			if (res[0] == 0)
			{
				noVarRows.add(x);
			}
			else if (countExpressed < minSamplesExpressed)
			{
				noVarRows.add(x);
			}
		}

		if (noVarRows.size() == 0)
		{
			System.out.println("There are no rows without variance");
			return;
		}

		MyMatrix remainder = new MyMatrix(	this.rows() - noVarRows.size(),
											this.cols());
		int nextRow = 0;
		int n = 0;
		int skip = noVarRows.get(n);
		n++;
		for (int x = 0; x < this.rows(); x++)
		{
			if (skip == x)
			{
				if (n < noVarRows.size())
					skip = noVarRows.get(n);
				n++;
				continue;
			}
			remainder.setRow(	nextRow,
								this.getRowHeaders()[x],
								this.getRowValues(x));
			nextRow++;
		}
		this.matrix = remainder.matrix;
		this.setRowHeaders(remainder.getRowHeaders());

		if (removedGenesFN != null)
		{
			JuhaPCA.PCA.log(" 5.1 Writing file from which genes without variance are removed");
			this.write(removedGenesFN);
		}
	}

	private int[] hasVariance(	int row,
								double minimum)
	{
		double firstValue = this.matrix.get(row,
											0);
		int[] results = new int[2];
		for (int c = 0; c < this.cols(); c++)
		{
			if (firstValue != this.matrix.get(	row,
												c))	;
			results[0] = 1;
			if (this.matrix.get(row,
								c) >= minimum)
				results[1]++;
		}
		return results;
	}

	public void correctForTotalReadCount(	double endCounts,
											double add)
	{
		//according to http://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29 (voom)
		//log 2 transformation is done after/outside this function though
		for (int c = 0; c < this.cols(); c++)
		{
			double total = 0;
			for (int r = 0; r < this.rows(); r++)
			{
				total += this.matrix.get(	r,
											c);
			}
			total++;
			for (int r = 0; r < this.rows(); r++)
			{
				this.matrix.set(r,
								c,
								(this.matrix.get(	r,
													c)
										+ add) / total * endCounts);//(geneExpression+0.5)/totalExpresion*10^6
			}
		}
	}

	public void log2Transform(double add)//add +1 before transforming
	{
		logTransform(	2,
						add);
	}

	public void pow(int i,
					double minus)
	{
		for (int x = 0; x < this.rows(); x++)
		{
			for (int y = 0; y < this.cols(); y++)
			{
				this.matrix.set(x,
								y,
								Math.pow(	10,
											this.matrix.get(x,
															y))
										- minus);
			}
		}
	}

	public void divideByCol(MatrixStruct denominators,
							int col)
	{
		for (int c = 0; c < this.cols(); c++)
		{
			for (int r = 0; r < this.rows(); r++)
			{
				this.matrix.set(r,
								c,
								this.matrix.get(r,
												c)
										/ denominators.matrix.get(	c,
																	col));
			}
		}
	}

	public void removeRow(String removeGene)
	{
		if (!this.getRowHash().containsKey(removeGene))
		{
			System.out.println("This matrix does not contain this gene in the rowheaders: " + removeGene);
			return;
		}
		else
			System.out.println("Removing: " + removeGene);
		MyMatrix result = new MyMatrix(	this.rows() - 1,
										this.cols());
		result.setColHeaders(this.getColHeaders());
		int outR = 0;
		for (int r = 0; r < this.rows(); r++)
		{
			if (this.getRowHeaders()[r].compareTo(removeGene) == 0)
			{
				continue;
			}
			result.setRow(	outR,
							this.getRowHeaders()[r],
							this.getRowValues(r));
			outR++;
		}
		this.rowNames = result.getRowHeaders();
		this.matrix = result.matrix;
	}

	public void putGenesOnRows()
	{
		if (this.getColHeaders()[0].contains("ENSG0") || this.getColHeaders()[0].contains("ENST0"))
		{
			System.out.println("Ensemble IDs are on columns and should be on rows, transposing");
			this.transpose();
		}
	}

	public void putGenesOnCols()
	{
		if (this.getRowHeaders()[0].contains("ENSG0") || this.getRowHeaders()[0].contains("ENST0"))
		{
			System.out.println("Ensemble IDs are on rows and should be on columns, transposing");
			this.transpose();
		}
	}

	public void removeRows(MatrixStruct removeGenes)
	{
		HashMap<String, Integer> toKeep = new HashMap<String, Integer>();
		int newPos = 0;
		for (int r = 0; r < this.rows(); r++)
		{
			if (!removeGenes.getRowHash().containsKey(this.getRowHeaders()[r]))
			{
				toKeep.put(	this.getRowHeaders()[r],
							newPos);
				newPos++;
			}
		}
		//remove the redundant IDs from the matrix
		this.keepIDs(toKeep);
	}

	public List<String> getRowHeadersAsList()
	{
		List<String> list = new ArrayList<String>();
		for (String rowName : this.getRowHeaders())
			list.add(rowName);
		return list;
	}

	public List<String> getColHeadersAsList()
	{
		List<String> list = new ArrayList<String>();
		for (String colName : this.getColHeaders())
			list.add(colName);
		return list;
	}

	public void divideColBy(int c,
							int denominator)
	{
		for (int r = 0; r < this.rows(); r++)
			this.matrix.set(r,
							c,
							this.matrix.get(r,
											c)
									/ denominator);
	}

	public void setColValues(	double[] values,
								int c)
	{
		for (int r = 0; r < this.rows(); r++)
			this.matrix.set(r,
							c,
							values[r]);
	}

	public Pair<Integer, Double> getRowLargest(int r)
	{
		double largest = 0;
		int index = 0;
		for (int c = 0; c < this.cols(); c++)
		{
			double value = this.matrix.get(	r,
											c);
			if (value > largest)
			{
				index = c;
				largest = value;
			}
		}
		Pair<Integer, Double> pair = new Pair<Integer, Double>(	index,
																largest);
		return pair;
	}

	public Pair<Integer, Double> getRowSmallest(int r)
	{
		double smallest = Double.POSITIVE_INFINITY;
		int index = 0;
		for (int c = 0; c < this.cols(); c++)
		{
			double value = this.matrix.get(	r,
											c);
			if (value < smallest)
			{
				index = c;
				smallest = value;
			}
		}
		Pair<Integer, Double> pair = new Pair<Integer, Double>(	index,
																smallest);
		return pair;
	}

	public MyMatrix mult(	MyMatrix B,
							MyMatrix C)
	{
		checkMultAdd(	B,
						C);
		for (int i = 0; i < this.rows(); ++i)
			for (int j = 0; j < C.cols(); ++j)
			{
				double dot = 0;
				for (int k = 0; k < this.cols(); ++k)
					dot += this.values[i][k] * B.values[k][j];
				C.values[i][j] += dot;
			}
		return C;
	}

	 public static MatrixStruct[] scores(MatrixStruct eigenVectors, MatrixStruct expression, String originalPath, boolean rotateBack, boolean correctForAverages) throws IOException 
	    {
	    	//putInCorrectOrientation(eigenVectors, expression, rotateBack);
	    	
	    	DenseMatrix originalMatrix = expression.matrix;
	    	String[] rowHeaders = eigenVectors.getRowHeaders();
	    	String[] colHeaders = eigenVectors.getColHeaders();
	    	DenseMatrix evMatrix = eigenVectors.matrix;
	    	int numRowsOri = expression.getRowHeaders().length;
	    	int numColsOri = expression.getColHeaders().length;
	    	int numRowsEV = eigenVectors.getRowHeaders().length;
	    	int numColsEV = eigenVectors.getColHeaders().length;
	    	
	    	log("Correcting for row averages");
	    	MatrixStruct[] scoreResults = new MatrixStruct[4];
	        log("Calculating principal component scores");
	        DenseMatrix scoreMatrix = new DenseMatrix(numRowsEV, numColsOri);
	        evMatrix.mult(originalMatrix, scoreMatrix);
	        log("Principal component scores calculated");
	        String[] pcHeaders = new String[numRowsEV];
	        for (int i = 0, len = pcHeaders.length; i < len; i++) {
	            pcHeaders[i] = "PC" + (i + 1);
	        }
	        scoreResults[0] = new MatrixStruct(scoreMatrix, pcHeaders,expression.getColHeaders());
	        if(rotateBack)//transpose
	        {
	        	scoreResults[0].transpose();
	        	scoreResults[0].setRowHeaders(expression.getColHeaders());//the rowheaders of the file that was rotated back become gene names again (which are on the columnHeaders of the eigenVector matrix)
	        }
	        log("Writing scores");
	        String FN = originalPath.toString().replace(".gz", "").replace(".txt", "") + ".scores.txt";
	        scoreResults[0].write(FN);

	        log("Calculating Cronbach's alpha for each component");
	        double[] alphas = cronbachsAlpha(evMatrix, scoreMatrix);
	        log("Cronbach's alphas calculated");

	        log("Writing Cronbach's alphas");
	        try (FileWriter fw = new FileWriter(originalPath.toString().replace(".gz", "").replace(".txt", "") + ".cronbachsAlpha.txt")) {
	            FileUtil.writeArray(fw, "cronbachsAlpha", alphas);
	        }
	        
	        scoreResults[2] = new MatrixStruct(originalMatrix, expression.getRowHeaders(), expression.getColHeaders());
	        scoreResults[3] = eigenVectors;
	        return scoreResults;
	    }
	
	protected void checkMultAdd(MyMatrix B,
								MyMatrix C)
	{
		if (rows() != C.rows())
			throw new IndexOutOfBoundsException("A.numRows != C.numRows (" + rows() + " != " + C.rows() + ")");
		if (cols() != B.rows())
			throw new IndexOutOfBoundsException("A.numColumns != B.numRows (" + cols() + " != " + B.rows() + ")");
		if (B.cols() != C.cols())
			throw new IndexOutOfBoundsException("B.numColumns != C.numColumns (" + B.rows() + " != " + C.cols() + ")");
	}

	public static double[] cronbachsAlpha(DenseMatrix evMatrix, DenseMatrix scoreMatrix) {

        int numComps = evMatrix.numRows();
        int len = evMatrix.numColumns();
        int lenScores = scoreMatrix.numColumns();

        double[] alphas = new double[numComps];
        for (int comp = 0; comp < numComps; comp++) {

            double evSquaredSum = 0;
            for (int i = 0; i < len; i++) {
                evSquaredSum += evMatrix.get(comp, i) * evMatrix.get(comp, i);
            }

            double scoreMean = 0;
            for (int i = 0; i < lenScores; i++) {
                scoreMean += scoreMatrix.get(comp, i) / lenScores;
            }

            double scoreVariance = 0;
            for (int i = 0; i < lenScores; i++) {
                double score = scoreMatrix.get(comp, i);
                scoreVariance += (score - scoreMean) * (score - scoreMean) / (lenScores - 1);
            }

            double alpha = (lenScores / (lenScores - 1d)) * (1d - (evSquaredSum / scoreVariance));
            alphas[comp] = alpha;
        }

        return alphas;
    }
	
	public MyMatrix toRatiosPerRow(boolean nanToZero)
	{
		MyMatrix ratios = new MyMatrix(	this.rows(),
										this.cols());
		ratios.rowNames = this.rowNames;
		ratios.colNames = this.colNames;
		ratios.firstField = this.firstField;
		
		for (int r = 0; r < rows(); r++)
		{
			double sum = 0;
			for (double val : this.values[r])
			{
				sum += val;
			}
			
			for (int c = 0; c < this.cols(); c++)
			{
				if (sum == 0)
					if (nanToZero)
					{
						ratios.values[r][c] = 0;
						continue;
					}
				ratios.values[r][c] = this.values[r][c] / sum;
			}
		}
		return ratios;
	}
	public MyMatrix toRatiosPerCol(boolean nanToZero)
	{
		MyMatrix ratios = new MyMatrix(	this.rows(),
										this.cols());
		ratios.rowNames = this.rowNames;
		ratios.colNames = this.colNames;
		ratios.firstField = this.firstField;
		
		for (int c = 0; c < this.cols(); c++)
		{
			double sum = 0;
			for (int r = 0; r < this.rows(); r++)
			{
				sum += this.values[r][c];
			}
			
			for (int r = 0; r < this.rows(); r++)
			{
				if (sum == 0)
					if (nanToZero)
					{
						ratios.values[r][c] = 0;
						continue;
					}
				ratios.values[r][c] = this.values[r][c] / sum;
			}
		}
		return ratios;
	}
	
    public static void log(String text) {

        String time = timeFormat.format(new Date());
        System.out.println(time + "\t" + text);
    }

	public MyMatrix mergeRows(MyMatrix matrixToAdd)
	{
		//keeps columns from the original matrix only
		this.getColHash();
		
		MyMatrix newMatrix = new MyMatrix(this.rows()+matrixToAdd.rows(), this.cols()+matrixToAdd.cols());
		
		newMatrix.colNames=this.colNames;
		for(int r=0; r< matrixToAdd.rows(); r++)
			newMatrix.rowNames[r]=this.rowNames[r];
		
		for(int r=0; r< this.rows(); r++)
			for(int c =0; c < this.cols(); c++)
				newMatrix.values[r][c]=this.values[r][c];


		//add the fields in the addmatrix that are also in this matrix
		
		//first rownames
		int newRow=this.rows();	
		for(int r=0; r< matrixToAdd.rows(); r++)
		{
			newMatrix.rowNames[newRow]=matrixToAdd.rowNames[r];
			newRow++;
		}
		
		for(int c =0; c < matrixToAdd.cols(); c++)
		{
			String colNameAddMatrix=matrixToAdd.colNames[c];
			
			if(!this.getColHash().containsKey(colNameAddMatrix))
				continue;
			
			int newCol=this.getColHash().get(colNameAddMatrix);
			newRow=this.rows();	
			for(int r=0; r< matrixToAdd.rows(); r++)
			{
				newMatrix.values[newRow][newCol]=this.values[r][c];
				newRow++;
				
			}
		}
		
		this.rowNames=newMatrix.rowNames;
		this.colNames=newMatrix.colNames;
		this.values=newMatrix.values;
		
		return newMatrix;
	}

	public String getRowInfoInOneLine(int row)
	{
		String line = this.rowNames[row];
		for(int c = 0; c< this.cols(); c++)
		{
			line+="\t"+this.colNames[c]+"="+this.values[row][c];
		}
		return line;
	}
}
