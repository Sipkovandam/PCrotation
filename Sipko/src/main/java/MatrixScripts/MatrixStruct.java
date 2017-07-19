package MatrixScripts;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.Hashtable;

import no.uib.cipr.matrix.DenseMatrix;
import umcg.genetica.containers.Pair;

import java.lang.Math;
import java.util.List;
import java.util.Random;
import java.util.Vector;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.stat.inference.WilcoxonSignedRankTest;

import JuhaPCA.FileUtil;
import JuhaPCA.PCA;
import org.apache.commons.math3.stat.inference.TTest;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.moment.Variance;

public class MatrixStruct
{
	//A java class that is limited to a Java array max size, but faster then Matrix

	public Hashtable<String, Integer> rowHash;
	public Hashtable<String, Integer> colHash;

	private String[] rowHeaders;
	private String[] colHeaders;
	public DenseMatrix matrix;

	public MatrixStruct(DenseMatrix matrix)
	{
		this.matrix = matrix;
	}

	public MatrixStruct(int nRows,
						int nCols)
	{
		this.setRowHeaders(makeHeaders(	"Row",
										nRows));
		this.setColHeaders(makeHeaders(	"Col",
										nCols));
		this.matrix = new DenseMatrix(	getRowHeaders().length,
										getColHeaders().length);
	}

	public MatrixStruct(DenseMatrix matrix,
						String[] rowHeaders,
						String[] colHeaders)
	{
		this.setRowHeaders(rowHeaders);
		this.setColHeaders(colHeaders);
		this.matrix = matrix;
	}

	public MatrixStruct(String[] rowNames,
						String[] colNames,
						double[][] expression)
	{
		this.setRowHeaders(rowNames);
		this.setColHeaders(colNames);
		this.setMatrix(expression);
	}

	public MatrixStruct(String fileName)
	{
		this(fileName, -1, -1);
	}

	public MatrixStruct(String fileName,
						int nRows,
						int nCols)
	{
		readFile(	fileName,
					nRows,
					nCols);

	}

	private String[] makeHeaders(	String prefix,
									int nCols)
	{
		String[] headers = new String[nCols];
		for (int h = 0; h < headers.length; h++)
			headers[h] = prefix + Integer.toString(h + 1);
		return headers;
	}

	private void setMatrix(double[][] expression)
	{
		matrix = new DenseMatrix(	expression.length,
									expression[0].length);
		for (int r = 0; r < expression.length; r++)
		{
			for (int c = 0; c < expression[0].length; c++)
			{
				matrix.set(	r,
							c,
							expression[r][c]);
			}
		}
	}

	void readFile(String fileName)
	{
		readFile(	fileName,
					-1,
					-1);
	}

	private void readFile(	String fileName,
							int nRows,
							int nCols)
	{
		Path path = Paths.get(fileName);
		try
		{
			setRowHeaders(FileUtil.readRowHeaders(	path,
													nRows));
			setColHeaders(FileUtil.readColumnHeaders(	path,
														nCols));
			nRows = this.getRowHeaders().length;
			nCols = this.getColHeaders().length;
			matrix = new DenseMatrix(	nRows,
										nCols);
			FileUtil.readMatrix(path,
								matrix,
								false,
								nRows,
								nCols);
		} catch (IOException e)
		{
			if (!new File(fileName).exists())
			{
				PCA.log("File does not exist: " + fileName);
				PCA.log("Exiting");
				e.printStackTrace();
			}
		}
	}

	public static Hashtable<String, Integer> makeHash(String[] headers)
	{
		if (headers == null)
			return null;
		Hashtable<String, Integer> temp = new Hashtable<String, Integer>();
		for (int x = 0; x < headers.length; x++)
		{
			if (headers[x] == null)
				continue;
			temp.put(	headers[x],
						x);
		}
		return temp;
	}

	public void transpose()
	{
		String[] temp = getRowHeaders();
		setRowHeaders(getColHeaders());
		setColHeaders(temp);

		DenseMatrix tempMat = new DenseMatrix(	matrix.numColumns(),
												matrix.numRows());

		for (int r = 0; r < matrix.numRows(); r++)
		{
			for (int c = 0; c < matrix.numColumns(); c++)
			{
				tempMat.set(c,
							r,
							matrix.get(	r,
										c));
			}

		}
		matrix = tempMat;
	}

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
			double average = rowAverages.matrix.get(rowAverages.getRowHash()
					.get(this.getRowHeaders()[r]),
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

	public void write(String fileName) throws IOException
	{
		write(fileName, false);
	}
	
	public void write(String fileName, boolean append) throws IOException
	{
		if (fileName.endsWith(".gz"))
		{
			GZIPOutputStream gzipOutputStream = new GZIPOutputStream(new FileOutputStream(fileName,append));

			try (BufferedWriter fw = new BufferedWriter(new OutputStreamWriter(gzipOutputStream)))
			{
				FileUtil.writeMatrix(	fw,
										matrix,
										getRowHeaders(),
										getColHeaders(),
										append);
			}
		}
		else
			try (FileWriter fw = new FileWriter(fileName,append))
			{
				FileUtil.writeMatrix(	fw,
										matrix,
										getRowHeaders(),
										getColHeaders(),
										append);
			}

		//if big matrix write smaller versions as well
		//		if(fileName.endsWith(".txt") && this.rows() > 1000 && this.cols()>1000)
		//		{
		//			//write first 10 rows
		//			try (FileWriter fw = new FileWriter(fileName.replace(".txt", "_10Rows.txt"))) {
		//				MatrixStruct subMatrix = this.getSubMatrix(10,-1);
		//	            FileUtil.writeMatrix(fw, subMatrix.matrix, subMatrix.getRowHeaders(), subMatrix.getColHeaders());
		//	        }
		//			//write first 10 columns
		//			try (FileWriter fw = new FileWriter(fileName.replace(".txt", "_10Cols.txt"))) {
		//				MatrixStruct subMatrix = this.getSubMatrix(-1,10);
		//	            FileUtil.writeMatrix(fw, subMatrix.matrix, subMatrix.getRowHeaders(), subMatrix.getColHeaders());
		//	        }
		//		}
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
			if (this.getRowHash()
					.get(PC) == null)
			{
				System.out.println("PC not found: " + PC);
				continue;
			}
			int r = this.getRowHash()
					.get(PC);
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

	public void setRowHeaders(String[] rowHeaders)
	{
		this.rowHeaders = rowHeaders;
		rowHash = makeHash(this.getRowHeaders());
	}

	public void setColHeaders(String[] colHeaders)
	{
		this.colHeaders = colHeaders;
		colHash = makeHash(this.getColHeaders());
	}

	public int rows()
	{
		return this.getRowHeaders().length;
	}

	public int cols()
	{
		return this.getColHeaders().length;
	}

	public String[] getRowHeaders()
	{
		return this.rowHeaders;
	}

	public String[] getColHeaders()
	{
		return this.colHeaders;
	}

	public MatrixStruct copy()
	{
		MatrixStruct copy = new MatrixStruct(	1,
												1);
		copy.setRowHeaders(deepStringCopy(this.getRowHeaders()));
		copy.setColHeaders(deepStringCopy(this.getColHeaders()));
		copy.matrix = new DenseMatrix(	this.matrix,
										true);//deepcopy
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
		this.rowHeaders[rowNumber] = rowName;
		this.rowHash = null;
		for (int c = 0; c < rowValues.length; c++)
		{
			this.matrix.set(rowNumber,
							c,
							rowValues[c]);
		}

	}
	public void setCol(	int colNumber,
						String colName,
						double[] values)
	{
		this.setColHeader(colNumber, colName);
		for (int r = 0; r < values.length; r++)
		{
			this.matrix.set(r,
			                colNumber,
							values[r]);
		}

	}

	public void keepRows1Matrix(MatrixStruct sample)
	{
		//keeps only the IDs in Sample
		Hashtable<String, Integer> toKeep = new Hashtable<String, Integer>();
		int newPos = 0;
		for (int r = 0; r < this.rows(); r++)
		{
			if (sample.getRowHash()
					.containsKey(this.rowHeaders[r]))
			{
				toKeep.put(	this.rowHeaders[r],
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
		Hashtable<String, Integer> toKeep = new Hashtable<String, Integer>();
		int newPos = 0;
		for (int r = 0; r < this.rows(); r++)
		{
			if (sample.getRowHash()
					.containsKey(this.rowHeaders[r]))
			{
				toKeep.put(	this.rowHeaders[r],
							newPos);
				newPos++;
			}
		}
		//remove the redundant IDs from the matrix
		this.keepIDs(toKeep);
		sample.keepIDs(toKeep);
	}

	public Hashtable<String, Integer> getRowHash()
	{
		if (this.rowHash == null)
			this.rowHash = makeHash(this.getRowHeaders());
		return this.rowHash;
	}

	public void keepIDs(Hashtable<String, Integer> toKeep)
	{
		MatrixStruct temp = new MatrixStruct(	toKeep.size(),
												this.cols());
		for (int r = 0; r < this.rows(); r++)
		{
			if (toKeep.get(this.rowHeaders[r]) == null)
				continue;
			int row = toKeep.get(this.rowHeaders[r]);
			//			System.out.println(this.rowHeaders[r] + " tokeepsize = " + toKeep.size());
			temp.setRow(row,
						this.rowHeaders[r],
						this.getRowValues(r));
		}
		this.setRowHeaders(temp.getRowHeaders());
		this.matrix = temp.matrix;
	}

	public MatrixStruct stDevRows()
	{
		MatrixStruct stDev = new MatrixStruct(	this.rows(),
												1);
		stDev.setColHeaders(new String[]
		{ "stDev" });
		stDev.setRowHeaders(this.rowHeaders);

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
		stDev.setColHeaders(new String[]
		{ "stDevCols" });
		stDev.setRowHeaders(this.colHeaders);

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

	public void divideBy(	MatrixStruct denoms,
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

	public MatrixStruct mergeColumns(MatrixStruct addition)
	{
		return mergeColumns(addition,
							false);
	}

	public MatrixStruct mergeColumns(	MatrixStruct addition,
										boolean keepAll)//adds columns to a matrix
	{
		int n = 0;
		if (keepAll == true)
			n = this.rows();
		else
			for (int r = 0; r < this.rows(); r++)
			{
				if (!addition.getRowHash()
						.containsKey(this.getRowHeaders()[r]))//if this row does not exist in the other file skip it
					continue;
				n++;
			}

		MatrixStruct result = new MatrixStruct(	n,
												this.cols() + addition.cols());
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

		//add the new rows
		int rOut = 0;
		for (int r = 0; r < this.rows(); r++)
		{
			int cOut = 0;
			double[] row = new double[this.cols() + addition.cols()];
			if (!keepAll && !addition.getRowHash()
					.containsKey(this.getRowHeaders()[r]))//if this row does not exist in the other file skip it
				continue;
			for (int c = 0; c < this.cols(); c++)
			{
				row[cOut] = this.matrix.get(r,
											c);
				cOut++;
			}
			if (addition.getRowHash()
					.get(this.getRowHeaders()[r]) != null)//if it exists in the added file
			{
				int rowToGet = addition.getRowHash()
						.get(this.getRowHeaders()[r]);
				for (int c = 0; c < addition.cols(); c++)
				{
					row[cOut] = addition.matrix.get(rowToGet,
													c);
					cOut++;
				}
			}
			result.setRow(	rOut,
							this.getRowHeaders()[r],
							row);
			rOut++;
		}
		System.out.println((this.rows() - result.rows()) + " rows removed (IDs not present in the added file)");
		return result;
	}

	public void setColHeader(	int c,
								String name)
	{
		this.colHeaders[c] = name;
		this.colHash = null;
	}

	public void setRowHeader(	int r,
								String name)
	{
		this.rowHeaders[r] = name;
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
	public void expressionToRank(MatrixStruct vector)
	{
		expressionToRank(	vector,
							0);
	}

	public void expressionToRank(	MatrixStruct vector,
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
				vector = new MatrixStruct(	this.rows(),
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

				//System.out.println("x =" + x + "Initial value: " + column[x][0] + "vNew value: " + val + " outputRow: " + (int) column[x][1] + " quantVectorRow: " + x + " c: " + y + " rowname = " + this.rowNames[(int) column[x][1]]);
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

	public MatrixStruct quantileNormVector()
	{
		MatrixStruct sortedCols = new MatrixStruct(	this.rows(),
													this.cols());
		sortedCols.setRowHeaders(this.getRowHeaders());
		sortedCols.setColHeaders(this.getColHeaders());
		for (int y = 0; y < this.cols(); y++)
		{
			double[] col = new double[this.rows()];
			for (int x = 0; x < this.rows(); x++)
			{
				col[x] = this.matrix.get(	x,
											y);
			}
			Arrays.sort(col);

			for (int x = 0; x < this.rows(); x++)
			{
				sortedCols.matrix.set(	x,
										y,
										col[x]);
			}
		}
		MatrixStruct qNormVector = getAveragesPerRow(sortedCols);
		return qNormVector;
	}

	public MatrixStruct getAveragesPerRow()
	{
		return getAveragesPerRow(this);
	}

	public MatrixStruct getAveragesPerRow(MatrixStruct data)
	{
		MatrixStruct qNormVector = new MatrixStruct(data.rows(),
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
		qNormVector.setColHeaders(new String[]
		{ "Averages" });
		return qNormVector;
	}

	public MatrixStruct getAveragesPerCol()
	{
		return getAveragesPerCol(this);
	}

	public MatrixStruct getAveragesPerCol(MatrixStruct data)
	{
		MatrixStruct temp = new MatrixStruct(	data.cols(),
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

	public double getAverageCol(int col)
	{
		double avg = 0;
		double denominator = rows();
		for (int r = 0; r < rows(); r++)
		{
			if (Double.isFinite(this.matrix.get(r,
												col)))
				avg += this.matrix.get(	r,
										col);
			else
				denominator--;
		}
		if (denominator == 0)
			denominator = 1;
		avg /= denominator;

		return avg;
	}

	public void adjustForAverageAllGenes(MatrixStruct rowAverages) //can be called with null and will use the for the averages of the columns of this matrix
	{
		for (int r = 0; r < rows(); r++)
		{
			adjustForAverageRow(rowAverages,
								r);
		}
	}

	public void adjustForAverageAllSamples(MatrixStruct columnAverages) //can be called with null and will use the for the averages of the columns of this MatrixStruct
	{
		for (int y = 0; y < cols(); y++)
		{
			adjustForAverageCol(columnAverages,
								y);
		}
	}

	public void adjustForAverageCol(MatrixStruct columnAverages,
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

		MatrixStruct remainder = new MatrixStruct(	this.rows() - noVarRows.size(),
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
												c))
				;
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

	public void logTransform(	int val,
								double add)//add +1 before transforming
	{
		double logVal = Math.log(val);
		for (int x = 0; x < this.rows(); x++)
		{
			for (int y = 0; y < this.cols(); y++)
			{
				this.matrix.set(x,
								y,
								Math.log(this.matrix.get(	x,
															y)
										+ add) / logVal);
			}
		}
	}

	public double getAverageRow(int r)
	{
		double sum = 0;
		int denom = this.cols();
		for (int c = 0; c < cols(); c++)
		{
			if (Double.isFinite(this.matrix.get(r,
												c)))
				sum += this.matrix.get(	r,
										c);
			else
				denom--;
		}

		double avg = sum / denom;
		return avg;
	}

	double sumRow(int r)
	{
		double sum = 0;
		for (int c = 0; c < cols(); c++)
		{
			//if(Double.isFinite(this.matrix.get(r,c)))
			sum += this.matrix.get(	r,
									c);
		}
		return sum;
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

	public void removeNoVariance(String removedGenesFN) throws IOException
	{
		//important that the probes/gene/transcripts are on the rows (rowNames)
		System.out.println("Rows before = " + this.rows());
		ArrayList<Integer> noVarRows = new ArrayList<Integer>();
		for (int x = 0; x < this.rows(); x++)//Identify all rows that have no variance
		{
			double var = new Variance().evaluate(this.getRowValues(x));
			if (var == 0)
				noVarRows.add(x);
		}

		if (noVarRows.size() == 0)
		{
			System.out.println("There are no rows without variance");
			return;
		}
		System.out.println("Removing:" + noVarRows.size() + " rows");

		//keep only the rows that have variance
		double[][] remainder = new double[this.rows() - noVarRows.size()][];
		String[] remainderRowNames = new String[this.rows() - noVarRows.size()];

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
			remainder[nextRow] = this.getRowValues(x);
			remainderRowNames[nextRow] = this.getRowHeaders()[x];
			nextRow++;
		}
		this.setMatrix(remainder);
		this.setRowHeaders(remainderRowNames);

		if (removedGenesFN != null)
		{
			JuhaPCA.PCA.log(" 5.1 Writing file from which genes without variance are removed");
			this.write(removedGenesFN);
		}
		System.out.println("Rows after" + this.rows() + " values=" + this.matrix.numRows());
	}

	public void removeRow(String removeGene)
	{
		if (!this.getRowHash()
				.containsKey(removeGene))
		{
			System.out.println("This matrix does not contain this gene in the rowheaders: " + removeGene);
			return;
		}
		else
			System.out.println("Removing: " + removeGene);
		MatrixStruct result = new MatrixStruct(	this.rows() - 1,
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
		this.rowHeaders = result.rowHeaders;
		this.matrix = result.matrix;
	}

	public void putGenesOnRows()
	{
		if (this.colHeaders[0].contains("ENSG0") || this.colHeaders[0].contains("ENST0"))
		{
			System.out.println("Ensemble IDs are on columns and should be on rows, transposing");
			this.transpose();
		}
	}

	public void putGenesOnCols()
	{
		if (this.rowHeaders[0].contains("ENSG0") || this.rowHeaders[0].contains("ENST0"))
		{
			System.out.println("Ensemble IDs are on rows and should be on columns, transposing");
			this.transpose();
		}
	}

	public void removeRows(MatrixStruct removeGenes)
	{
		Hashtable<String, Integer> toKeep = new Hashtable<String, Integer>();
		int newPos = 0;
		for (int r = 0; r < this.rows(); r++)
		{
			if (!removeGenes.getRowHash()
					.containsKey(this.rowHeaders[r]))
			{
				toKeep.put(	this.rowHeaders[r],
							newPos);
				newPos++;
			}
		}
		//remove the redundant IDs from the matrix
		this.keepIDs(toKeep);
	}

	public MatrixStruct getRow(int row) // deep copy 
	{
		MatrixStruct matRow = new MatrixStruct(	this.cols(),
												1);
		matRow.setRowHeaders(new String[this.cols()]);
		matRow.setColHeaders(new String[]
		{ this.getRowHeaders()[row] });
		for (int c = 0; c < this.cols(); c++)
		{
			matRow.matrix.set(	c,
								0,
								this.matrix.get(row,
												c));
			matRow.setRowHeader(c,
								this.getColHeaders()[c]);
		}
		matRow.rowHash = makeHash(matRow.getRowHeaders());
		matRow.colHash = makeHash(matRow.getColHeaders());
		return matRow;
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

//	public void putGenesOnCorrectAxis(boolean isGenesOnRows)
//	{
//		if (isGenesOnRows)
//			putGenesOnRows();
//		else
//			putGenesOnCols();
//	}
}
