package PCA;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.LineNumberReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.Hashtable;
import java.nio.file.Path;

public class Matrix 
{
	public String firstField;
	public String[] colNames;
	public String[] rowNames;
	public double[][] values;
	public Matrix avgCols = null;
	public boolean verbose = false;
	
	
	public Matrix(int x, int y)
	{
		rowNames = new String[x];
		colNames = new String[y];
		values = new double[x][y];
		
		numberNames(rowNames, "Row");
		numberNames(colNames, "Col");		
	}
	
	public Matrix(int x, int y, String[] colNames, String[] rowNames)
	{
		this.rowNames = rowNames;
		this.colNames = colNames;
		values = new double[x][y];
	}
	
	
	public Matrix(String fileName)
	{
		readFile(fileName,true,true);
	}
	public Matrix(Path fileName)
	{
		readFile(fileName.toString(),true,true);
	}
	public Matrix(String fileName,int x, int y)//read file with a maximum of xRows and yCols
	{
		readFile(fileName,true,true, x, y);
	}
	
	public Matrix(String fileName, boolean hasRowColNames)
	{
		readFile(fileName,false,false);
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

	public void log2Transform()//add +1 before transforming
	{
		logTransform(2);
	}
	public void logTransform(int val)//add +1 before transforming
	{
		double logVal = Math.log(val);
		for(int x = 0; x < this.rowNames.length; x++)
		{
			for(int y = 0; y < this.colNames.length; y++)
			{
				//System.out.println(this.values[x][y]+ " " +  Math.log(this.values[x][y]) + " " + logVal);
				//if(values[x][y] ==0)
				//	continue;
				this.values[x][y] = Math.log(this.values[x][y]+1)/logVal;
				/*if(x < 5 && y < 5)
				{
					System.out.println(this.values[x][y] + " " + (Math.log(this.values[x][y])/logVal) + " " + this.values[x][y]);
				}*/
			}
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
			BufferedReader reader = new BufferedReader(new FileReader(new File(fileName)));
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
				values = new double[nRows-1][nCols-minusCol];
			else if (hasColNames && !hasRowNames)
				values = new double[nRows][nCols-minusCol];
			else if (!hasColNames && hasRowNames)
				values = new double[nRows-1][nCols];
			else
				values = new double[nRows][nCols];
			
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
				reader = new BufferedReader(new FileReader(new File(fileName)));
			}
			if(hasRowNames)
				rowNames = new String[nRows-1];
			else
				rowNames = new String[nRows];
			
			//create the array for the rowNames & values
			//System.out.println("rowlength = "+ rowNames.length);
			
			
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
						if(hasRowNames)
							values[x][y]= Double.parseDouble(eles[y+1]);
						else
							values[x][y]= Double.parseDouble(eles[y]);
					}catch (Exception e1)
					{
						System.out.println("FileName =  " + fileName);
						System.out.println("cols =  " + values[x].length+" y = " +  y + "cols.length  = " + colNames.length);
						System.out.println("exception =  " + e1);
						System.out.println("x =  " + x + " y = " + y + " cols[y] = " + colNames[y]);
						values[x][y] = 0;
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
		} catch (Exception e) {e.printStackTrace();System.out.println("FileName =" + fileName);}

		if(verbose)
			System.out.println("Finished reading file");
		
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
	
	
	
	public void keepRows(Matrix m2)//Changes both matrixes if necessary, Makes sure same IDs are on same line in both files
	{
		//find all the IDs shared by both files
		Hashtable<String,Integer> allIDs = rowNamesToHash();
		Hashtable<String,Integer> toKeep = new Hashtable<String,Integer>();
		int newPos = 0;
		for(int r = 0; r < m2.rowNames.length;r++)
		{
			if(allIDs.containsKey(m2.rowNames[r]))
			{
				toKeep.put(m2.rowNames[r], newPos);
				newPos++;
			}
		}
		keepIDsInHash(toKeep);
		m2.keepIDsInHash(toKeep);
	}
	
	private void keepIDsInHash(Hashtable<String, Integer> toKeep) 
	{
		int total = toKeep.size();
		
		double[][] vals = new double[total][];
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
			{
					writer = new BufferedWriter(new FileWriter(new File(fileName)));	
			}
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
				if(x%(100000000/maxY)==0 && x>0)//print time 
				{
					
					//double percentage = (((double)x)/((double)maxX));
					//System.out.println(df.format(percentage*100) + " % of the file has been saved");
					//timer.print(percentage);
				}
			}
			if(writer != null)
				writer.close();
		} catch (IOException e) {System.out.println("Oh noes! An error in your print function.");e.printStackTrace();}
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
			FileReader fileReader = new FileReader(new File(fileName));
			LineNumberReader lnr = new LineNumberReader(fileReader);
			lnr.skip(Long.MAX_VALUE);
			nLines = lnr.getLineNumber() + 1; //Add 1 because line index starts at 0
			fileReader.close();
		} catch (FileNotFoundException e1) {e1.printStackTrace();} catch (IOException e) {e.printStackTrace();}
		//System.out.println("lines = " + nLines);		
		return nLines;
	}
	public void removeNoVariance()
	{
		//important that he probes/gene/transcripts are on the X-axis (rowNames)
		ArrayList<Integer> noVarRows = new ArrayList<Integer>();
		for(int x = 0; x < this.rowNames.length;x++)//Identify all rows that have no variance
		{
			boolean hasVariance = hasVariance(x);
			if(!hasVariance)
			{
				noVarRows.add(x);
			}
		}
		
		if(noVarRows.size() == 0)
		{
			System.out.println("There are no rows without variance");
			return;
		}
		
		double[][] remainder = new double[this.rowNames.length-noVarRows.size()][];
		String[] remainderRowNames = new String[this.rowNames.length-noVarRows.size()];
		
		int nextRow = 0;
		int n = 0;
		int skip = noVarRows.get(n);n++;
		for(int x = 0; x < this.rowNames.length; x++)
		{
			if(skip == x)
			{
				if(n < noVarRows.size())
					skip = noVarRows.get(n);
				n++;
				continue;
			}
			remainder[nextRow]=this.values[x];
			remainderRowNames[nextRow] = this.rowNames[x]; 
			nextRow++;
		}
		this.values = remainder;
		this.rowNames = remainderRowNames;
	}
	
	
	private boolean hasVariance(int row) 
	{
		double firstValue = this.values[row][0];
		for(int y = 0; y < this.colNames.length; y++)
		{
			if(firstValue != this.values[row][y])
				return true;
		}
		return false;
	}

	public Matrix calcAvgCols()
	{
		Matrix temp = new Matrix(this.colNames.length, 1);
		temp.rowNames=this.colNames;
		temp.colNames[0] = "Average";
		for(int y = 0; y < colNames.length; y++)
		{
			temp.values[y][0] = getAverageCol(y);
		}
		return temp;
	}
	public Matrix calcAvgRows()
	{
		Matrix temp = new Matrix(this.rowNames.length, 1);
		temp.rowNames=this.rowNames;
		temp.colNames[0] = "Average";
		for(int r = 0; r < rowNames.length;r++)
		{
			temp.values[r][0] = getAverageRow(r);
		}
		return temp;
	}
	public double getAverageRow(int r)
	{
		double avg = 0;
		for(int c = 0; c < colNames.length; c++)
		{
			avg+=values[r][c];
		}
		avg /= colNames.length;
		return avg;
	}
	public double getAverageCol(int col)
	{
		double avg = 0;
		for(int x = 0; x < rowNames.length; x++)
		{
			avg+=values[x][col];
		}
		avg /= rowNames.length;
		
		return avg;
	}
	
	public void adjustForAverageCol(Matrix columnAverages, int col)
	{
		double average = 0;
		if(columnAverages == null)
			average = getAverageCol(col);
		else
			average = columnAverages.values[col][0];
		for(int x = 0; x < rowNames.length; x++)
		{
			//System.out.println("x =" + x + " col =" + col + " values[x][col]=" + values[x][col] + " average=" + average);
			values[x][col]-= average;
		}
	}
	public void adjustForAverageRow(Matrix rowAverages, int row)
	{
		double average = 0;
		Hashtable<String, Integer> averagesRowHash = null;
		if(rowAverages != null)
		{
			averagesRowHash = rowAverages.rowNamesToHash();
//			System.out.println("this.rowNames[row]"+this.rowNames[row]);
//			System.out.println("averagesRowHash.get(this.rowNames[row]) = "+averagesRowHash.get(this.rowNames[row]));
//			System.out.println("averagesRowHash size = "+averagesRowHash.size());
//			System.out.println("averagesRowHash SRR1028344 = "+averagesRowHash.containsKey("SRR1028344"));
//			System.out.println("averagesRowHash element ERR315329 = "+averagesRowHash.containsKey("ERR315329"));
//			System.out.println("averagesRowHash element = "+averagesRowHash.keys().nextElement());
//			Enumeration<String> key = averagesRowHash.keys();key.nextElement();
//			System.out.println("averagesRowHash element = "+key.nextElement());
//			System.out.println("rowAverages = "+rowAverages);
			
			average = rowAverages.values[averagesRowHash.get(this.rowNames[row])][0];
		}
		else
			average = getAverageRow(row);
		for(int c = 0; c < colNames.length; c++)
		{
			values[row][c]-= average;
		}
	}

	public void adjustForAverageAllrows(Matrix rowAverages) //can be called with null and will use the for the averages of the columns of this matrix
	{
		for(int r =0; r < rowNames.length; r++)
		{
			adjustForAverageRow(rowAverages, r);
		}
	}
	public void adjustForAverageAllCols(Matrix columnAverages) //can be called with null and will use the for the averages of the columns of this matrix
	{
		for(int y =0; y < colNames.length; y++)
		{
			adjustForAverageCol(columnAverages, y);
		}
	}
	public Matrix quantileNormVector()
	{
		Matrix sortedCols = new Matrix(this.rowNames.length, this.colNames.length);
		sortedCols.rowNames=rowNames;
		sortedCols.colNames=colNames;
		for(int y = 0; y < this.colNames.length; y++)
		{
			double[] col = new double[this.rowNames.length];
			for(int x = 0; x < this.rowNames.length; x++)
			{
				col[x] = this.values[x][y];
			}
			Arrays.sort(col);
			
			for(int x = 0; x < this.rowNames.length; x++)
			{
				sortedCols.values[x][y] = col[x];
			}
		}
		Matrix qNormVector = averagesPerRow(sortedCols);
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

	public Matrix quantileNormVectorCenteredLog2(String saveNameQnormVector)//Calculated the averages for each rank
	{
		Matrix qNormVector = quantileNormVector();
		qNormVector.log2Transform();
		qNormVector.adjustForAverageAllCols(null);
		qNormVector.write(saveNameQnormVector.replace(".txt", "_TransformedCentered.txt"));
		System.out.println("Quantile normalization vector saved to " + saveNameQnormVector);
		qNormVector.rowNames = this.rowNames;//just need these to determine which genes to include in quantile normalization.
		qNormVector.colNames[0] = "centered log2 transformed quantile vector. RowNames have nothing to do with the quant vector!";
		return qNormVector;
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

	public Matrix averagesPerRow(Matrix sortedCols) 
	{
		Matrix qNormVector = new Matrix(sortedCols.rowNames.length, 1);
		for(int x = 0; x < sortedCols.rowNames.length; x++)
		{
			double average = 0;
			for(int y = 0; y < sortedCols.colNames.length; y++)
			{
				average += sortedCols.values[x][y];
			}
			average /= sortedCols.colNames.length;
			
			qNormVector.values[x][0] = average;
		}
		qNormVector.rowNames = sortedCols.rowNames;
		qNormVector.colNames = new String[]{"Averages"};
		return qNormVector;
	}

	public void quantileNormAdjust(Matrix qNormVector)//this comes after calculating the determining the means for each rank
	{
		if(qNormVector.rowNames.length!= this.rowNames.length)
		{
			System.out.println("Dimensions incorrect, transposing");
			this.transpose();
			if(qNormVector.rowNames.length!= this.rowNames.length)
			{
				System.out.println("Dimensions still incorrect incorrect: \n"
						+ " this.rowNames.length = " + this.rowNames.length + " this.colNames.length " + this.rowNames.length +"/n"
						+ " qNormVector.rowNames.length = " + qNormVector.rowNames.length);
			}
		}
		
		double[][] column= new double[this.rowNames.length][2];
		for(int y = 0; y < this.colNames.length; y++)//sort each column
		{
			for(int x = 0; x < this.rowNames.length; x++)
			{
				column[x][0] = this.values[x][y];
				column[x][1] = x;
				//System.out.println(this.rowNames[x] +" "+ x + " val:" +this.values[x][y]);
			}
			
			Arrays.sort(column, new Comparator<double[]>()
			{
				public int compare(double[] s1, double[] s2)
				{
					if (s1[0] > s2[0])
						return 1;	// tells Arrays.sort() that s1 comes after s2
					else if (s1[0] < s2[0])
						return -1;   // tells Arrays.sort() that s1 comes before s2
					else 
					{
						return 0;
					}
				}
			});			
			
			int qNormVectorIndAdj = 0;//to adjust for values that are equal to another value in the same column
			int pos = 0;
			double val = 0;
			for(int x = 0; x < column.length; x++)
			{
				//far to complicated if else to achieve something simple, but its late :S	
				if(x>0 && column[x][0] != column[x-1][0] && (x+1<column.length && column[x][0] != column[x+1][0]))//if this value is not equal to next or previous value
				{
					qNormVectorIndAdj=0;
					pos = x;		
					val = qNormVector.values[pos][0];
				}
				else
				{
					if (x==0 || (x>0 && column[x][0] != column[x-1][0]))//if this value is different from the last one recalculate the number to use (or if it is the first value)
					{
						qNormVectorIndAdj=0;
						int extra = 1;
						while(x+extra<column.length && column[x][0] == column[x+extra][0])
						{
							qNormVectorIndAdj++;
							extra++;
						}
						pos = x+(qNormVectorIndAdj/2);
					}

					if(qNormVectorIndAdj % 2 == 0 || pos == qNormVector.values.length-1)
						val = qNormVector.values[pos][0];
					else
					{
						val = (qNormVector.values[pos][0]+qNormVector.values[pos+1][0])/2;
					}
				}

				//System.out.println("x =" + x + "Initial value: " + column[x][0] + "vNew value: " + val + " outputRow: " + (int) column[x][1] + " quantVectorRow: " + x + " c: " + y + " rowname = " + this.rowNames[(int) column[x][1]]);
				this.values[(int) column[x][1]][y] = val;		
			}
		}
		
	} 
	
	public void transform(Matrix mat)//transforms the values of this matrix to the new PC space. not to be confused with transpose ;)
	{
		//If the expression matrix has the samples on the Y-axis and the genes on the X-axis
		//If the vectors of the different PCs are in columns (where each factor is the PC for each gene)
		//multiply each column of the expression matrix with each row of the PC matrix
		Matrix temp = new Matrix(mat.rowNames.length, this.colNames.length);//rowNames should be the samples
		temp.colNames=this.colNames;
		//give the rowNames PC1, PC2, PC3, etc..
		for(int x = 0; x < temp.rowNames.length; x++)
		{
			temp.rowNames[x] = "PC" + Integer.toString(x+1);
		}
		Timer timer = new Timer();
		//Calculate the new cooridinate for 1 sample each iteration
		for(int y = 0; y < this.colNames.length; y++)
		{
			if(y % 1000 == 0)//(this.colNames.length/10000)
			{
				double percentage = ((double)(y)/((double)(this.colNames.length)));
				System.out.println( ((int)(percentage*100)) + "% completed");
				timer.print(percentage);
			}
			for(int x = 0; x < mat.rowNames.length; x++)
			{
				//System.out.println("x = " + x + "/" + mat.rowNames.length);
				temp.values[x][y] = multiplyVectors(mat, x, y);
			}
		}
		
		this.values=temp.values;
		this.rowNames = temp.rowNames;
		
		//resulting file (temp), has PCs on RowNames, samples as colNames
	}
	private double multiplyVectors(Matrix mat, int col, int y2) //multiply a column of the expression data with a row of the PC matrix (mat= PC matrix)
	{
		double result = 0;
		double[] product = new double[mat.colNames.length];
		for(int y = 0; y < mat.colNames.length; y++)
		{
			//System.out.println("col = " + col + " y = " + y+ " this.values[y][col] " + this.values[y][y2] + " y2 = " + y2);
			product[y]=this.values[y][y2]*mat.values[col][y];//x in values is the columns	
		}
		double sum = sumArray(product);
		result = sum;
		return result;
	}

	private double sumArray(double[] numbers) 
	{
		double sum = 0;
		for(int x = 0; x < numbers.length; x++)
		{
			sum+=numbers[x];
		}
		return sum;
	}
	public void transpose()
	{
		double[][] temp = new double[colNames.length][rowNames.length];
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
					
						writer = new BufferedWriter(new FileWriter(new File(fileName)));
					
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
			}
		}catch (IOException e)
		{
			e.printStackTrace();
		}	
	}

	public Matrix getCol(int i) 
	{
		Matrix temp = new Matrix(this.rowNames.length,1);
		temp.rowNames = this.rowNames;
		temp.colNames = new String[]{this.colNames[i]};
		for(int r = 0; r < this.rowNames.length; r++)
		{
			temp.values[r] = this.values[r];
		}
		// TODO Auto-generated method stub
		return temp;
	}

	public Matrix getRow(int i) 
	{
		Matrix temp = new Matrix(1,this.colNames.length);
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
}
