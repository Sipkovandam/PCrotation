package JuhaPCA;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;
import no.uib.cipr.matrix.DenseMatrix;

/**
 *
 * @author juha
 */
public class FileUtil {

    private static int rowsInserted = 0;

    public static int countLines(Path path) throws IOException {

        int numLines;
        if (path.toString().endsWith(".gz")) {
            numLines = (int) GZipFiles.lines(path).count();
        } else {
            numLines = (int) Files.lines(path).count();
        }
        return numLines;
    }
    public static String[] readRowHeaders(Path path) throws IOException {
    	return  readRowHeaders(path, -1);
    }
    public static String[] readRowHeaders(Path path, int limit) throws IOException {
    	if(limit < 0)limit = countLines(path);
        List<String> headers = new ArrayList<>();
        if (path.toString().endsWith(".gz")) {
            GZipFiles.lines(path)
                    .skip(1) // header
                    .limit(limit)
                    .map(line -> line.split("\\t"))
                    .forEach(row -> headers.add(row[0]));
        } else {
            Files.lines(path)
                    .skip(1) // header
                    .limit(limit)
                    .map(line -> line.split("\\t"))
                    .forEach(row -> headers.add(row[0]));
        }
        return headers.toArray(new String[0]);
    }
    public static String[] readColumnHeaders(Path path) throws IOException
    {
    	return readColumnHeaders(path, -1);
    }
    
    public static String[] readColumnHeaders(Path path, int nCols) throws IOException {

        String[] headers;
        if (path.toString().endsWith(".gz")) {
            headers = GZipFiles.lines(path).findFirst().get().split("\\t");
        } else {
            headers = Files.lines(path).findFirst().get().split("\\t");//if you use .firstFirst.toString it will add "Optional[" and "]"
        }
        if(nCols < 1) nCols = headers.length;
        if(nCols > headers.length) nCols = headers.length;
        headers = Arrays.copyOfRange(headers, 1, nCols);
        return headers;
    }

    public static void readMatrix(Path path, DenseMatrix matrix) throws IOException {

        FileUtil.readMatrix(path, matrix, false);
    }
    public static void readMatrix(Path path, DenseMatrix matrix, boolean transpose) throws IOException {

        FileUtil.readMatrix(path, matrix, false, -1,-1);
    }

    public static void readMatrix(Path path, DenseMatrix matrix, boolean transpose, int nRows, int nCols) throws IOException {
        rowsInserted = 0;
        //System.out.println(path.toString());
        if(nRows < 0)nRows = countLines(path);
        if (path.toString().endsWith(".gz")) {
            GZipFiles.lines(path)
                    .skip(1) // header
                    .limit(nRows)
                    .map(line -> line.split("\\t"))
                    .forEach(row -> insert(row, matrix, -1, transpose,nCols));
        } else {
            Files.lines(path)
                    .skip(1) // header
                    .limit(nRows)
                    .map(line -> line.split("\\t"))
                    .forEach(row -> insert(row, matrix, -1, transpose,nCols));
        }
    }

    public static void readArray(Path path, DenseMatrix matrix, int column) throws IOException {
        readArray(path, matrix, column, false);
    }
    
    public static void readArray(Path path, DenseMatrix matrix, int column, boolean transpose) throws IOException {

        rowsInserted = 0;
        if (path.toString().endsWith(".gz")) {
            GZipFiles.lines(path)
                    .skip(1) // header
                    .map(line -> line.split("\\t"))
                    .forEach(row -> insert(row, matrix, column, transpose));
        } else {
            Files.lines(path)
                    .skip(1) // header
                    .map(line -> line.split("\\t"))
                    .forEach(row -> insert(row, matrix, column, transpose));
        }
    }
    private static void insert(String[] row, DenseMatrix matrix, int column, boolean transpose)
    {
    	insert(row, matrix, column, transpose, -1);
    }
    private static void insert(String[] row, DenseMatrix matrix, int column, boolean transpose, int nCols) {

        if (column > -1) {//reading only 1 column in the file
            if (transpose) {
                matrix.set(0, rowsInserted, Double.parseDouble(row[column]));
            } else {
                matrix.set(rowsInserted, 0, Double.parseDouble(row[column]));
            }
        } else {
        	if(nCols < 1) nCols = row.length;
        	else
        		nCols++;
            for (int i = 1, len = nCols; i < len; i++) {
                if (transpose) {
                    matrix.set(i - 1, rowsInserted, Double.parseDouble(row[i]));
                } else {
                    matrix.set(rowsInserted, i - 1, Double.parseDouble(row[i]));
                }
            }
        }
        ++rowsInserted;
    }

    public static void writeMatrix(FileWriter fw, DenseMatrix matrix, String[] rowHeaders, String[] colHeaders) throws IOException {
        writeMatrix(fw, matrix, rowHeaders, colHeaders, false);
    }
    public static void writeMatrix(BufferedWriter fw, DenseMatrix matrix, String[] rowHeaders, String[] colHeaders) throws IOException {
        writeMatrix(fw, matrix, rowHeaders, colHeaders, false);
    }
    public static void writeMatrix(FileWriter fw, DenseMatrix matrix, String[] rowHeaders, String[] colHeaders, boolean transpose) throws IOException {
    	writeMatrix(new BufferedWriter(fw), matrix, rowHeaders, colHeaders, transpose);
    }

    public static void writeMatrix(BufferedWriter fw, DenseMatrix matrix, String[] rowHeaders, String[] colHeaders, boolean transpose) throws IOException {
        String date = PCA.dateFormat.format(new Date());

        // write date and column headers to first line
        fw.write(date);
        if (transpose) {
            for (int c = 0, cols = matrix.numRows(); c < cols; c++) {
                fw.write("\t" + colHeaders[c]);
            }
        } else {
            for (int c = 0, cols = matrix.numColumns(); c < cols; c++) {
                fw.write("\t" + colHeaders[c]);
            }
        }
        fw.write("\n");

        // write data with row headers
        if (transpose) {
            for (int c = 0, cols = matrix.numColumns(); c < cols; c++) {
                fw.write(rowHeaders[c]);
                for (int r = 0, rows = matrix.numRows(); r < rows; r++) {
                    fw.write("\t" + matrix.get(r, c));
                }
                fw.write("\n");
            }
        } else {
            for (int r = 0, rows = matrix.numRows(); r < rows; r++) {
                fw.write(rowHeaders[r]);
                for (int c = 0, cols = matrix.numColumns(); c < cols; c++) {
                    fw.write("\t" + matrix.get(r, c));
                }
                fw.write("\n");
            }
        }
        fw.close();
    }

    public static void writeArray(FileWriter fw, String name, double[] array) throws IOException {

        fw.write(name + "\n");
        for (int i = 0, len = array.length; i < len; i++) {
            fw.write(array[i] + "\n");
        }
    }

}
