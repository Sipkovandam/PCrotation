package JuhaPCA;

import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.text.Format;
import java.text.SimpleDateFormat;
import java.util.Date;

import no.uib.cipr.matrix.DenseMatrix;


//code taken from Juha Karjalainen: https://github.com/juhis/PCA/blob/master/src/pca/PCA.java
public class CronbachCalculator
{
	static final Format timeFormat = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
	
	public static void cronbach(Path originalPath, Path evPath, Path scorePath, boolean isEVTransposed) throws IOException {

        if (isEVTransposed) {
            log("Reading eigenvectors (transposed)");
        } else {
            log("Reading eigenvectors");
        }
        int numRows = FileUtil.readColumnHeaders(evPath).length;
        log("Matrix size " + numRows + " x " + numRows);
        DenseMatrix evMatrix = new DenseMatrix(numRows, numRows);
        FileUtil.readMatrix(evPath, evMatrix, isEVTransposed);
        log("Eigenvectors read");

        log("Reading scores");
        String[] colHeaders = FileUtil.readColumnHeaders(scorePath);
        int numCols = colHeaders.length;
        DenseMatrix scoreMatrix = new DenseMatrix(numRows, numCols);
        FileUtil.readMatrix(scorePath, scoreMatrix);
        log("Scores read");

        log("Reading original data");
        String[] headers = FileUtil.readColumnHeaders(originalPath);
        numCols = headers.length;
        DenseMatrix originalMatrix;
        if (numCols == numRows) { // transpose data
            headers = FileUtil.readRowHeaders(originalPath);
            numCols = headers.length;
            log("Matrix size (transposed) " + numRows + " x " + numCols);
            originalMatrix = new DenseMatrix(numRows, numCols);
            FileUtil.readMatrix(originalPath, originalMatrix, true);
        } else {
            log("Matrix size " + numRows + " x " + numCols);
            originalMatrix = new DenseMatrix(numRows, numCols);
            FileUtil.readMatrix(originalPath, originalMatrix);
        }
        log("Original data read");

        log("Calculating Cronbach's alpha for each component");
        double[] alphas = cronbachsAlpha(originalMatrix, evMatrix, scoreMatrix);
        log("Cronbach's alphas calculated");

        log("Writing Cronbach's alphas");
        try (FileWriter fw = new FileWriter(scorePath.toString().replace(".gz", "").replace(".txt", "") + ".cronbachsAlpha.txt")) {
            FileUtil.writeArray(fw, "cronbachsAlpha", alphas);
        }
    }
	
	/**
    *
    * Calculate Cronbach's alpha for each principal component.
    *
    * @param evMatrix Eigenvector matrix, each row is an eigenvector
    * @param scoreMatrix Principal component score matrix, each row is a
    * component
    * @return An array of Cronbach's alpha values
    */
   private static double[] cronbachsAlpha(DenseMatrix originalMatrix, DenseMatrix evMatrix, DenseMatrix scoreMatrix) {

       int numComps = evMatrix.numRows();
       int lenItems = evMatrix.numColumns();
       int lenScores = scoreMatrix.numColumns();

       double[] variances = new double[lenItems];
       for (int i = 0; i < lenItems; i++) {
           variances[i] = getRowVariance(originalMatrix, i);
       }

       double[] alphas = new double[numComps];
       for (int comp = 0; comp < numComps; comp++) {

           double sumVariance = 0;
           for (int i = 0; i < lenItems; i++) {
               sumVariance += evMatrix.get(i, comp) * evMatrix.get(i, comp) * variances[i];
           }

           double scoreMean = getRowMean(scoreMatrix, comp);
           double scoreVariance = getRowVariance(scoreMatrix, scoreMean, comp);

           double alpha = (lenItems / (lenItems - 1d)) * (1d - (sumVariance / scoreVariance));
           alphas[comp] = alpha;
       }

       return alphas;
   }
   
   /**
   *
   * Calculates mean of a given row.
   *
   * @param matrix Data matrix
   * @param row Row index
   * @return Arithmetic mean of the row
   */
  private static double getRowMean(DenseMatrix matrix, int row) {

      double mean = 0;
      for (int c = 0, numCols = matrix.numColumns(); c < numCols; c++) {
          mean += matrix.get(row, c) / numCols;
      }
      return mean;
  }
  
   /**
   *
   * Calculates variance for a given row.
   *
   * @param matrix Data matrix
   * @param mean Mean of the row
   * @param row Row index
   * @return Sample variance of the row
   */
  private static double getRowVariance(DenseMatrix matrix, double mean, int row) {

      double variance = 0;
      for (int c = 0, numCols = matrix.numColumns(); c < numCols; c++) {
          variance += (matrix.get(row, c) - mean) * (matrix.get(row, c) - mean) / (numCols - 1);
      }
      return variance;
  }

  /**
   *
   * Calculates variance for a given row.
   *
   * @param matrix Data matrix
   * @param row Row index
   * @return Sample variance of the row
   */
  private static double getRowVariance(DenseMatrix originalMatrix, int row) {

      double mean = getRowMean(originalMatrix, row);
      return getRowVariance(originalMatrix, mean, row);
  }
   
    /**
    *
    * Print given text to stdout, preceded by a MySQL-style timestamp and a
    * tab.
    *
    * @param text Text to log
    */
   private static void log(String text) {

       String time = timeFormat.format(new Date());
       System.out.println(time + "\t" + text);
   }
}
