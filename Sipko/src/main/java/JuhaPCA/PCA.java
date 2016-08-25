package JuhaPCA;

import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.Format;
import java.text.SimpleDateFormat;
import java.util.Date;

import PCA.MatrixStruct;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.EVD;
import no.uib.cipr.matrix.NotConvergedException;

public class PCA {

    static final Format dateFormat = new SimpleDateFormat("yyyy-MM-dd");
    static final Format timeFormat = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

    /**
     *
     * Run one of the defined matrix operations:
     *
     * center, scale, transpose, covariance, correlation, evd, scores, transform
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        if (args.length < 1) {
            printUsage();
            System.exit(1);
        }

        if ("evd".equals(args[0]) && args.length != 2) {
            printUsage();
            System.exit(1);
        }

        if ("scores".equals(args[0]) && (args.length != 3 && args.length != 4)) {
            printUsage();
            System.exit(1);
        }

        if ("transpose".equals(args[0]) && args.length != 2) {
            printUsage();
            System.exit(1);
        }

        if ("center".equals(args[0]) && args.length != 2) {
            printUsage();
            System.exit(1);
        }

        if ("scale".equals(args[0]) && args.length != 2) {
            printUsage();
            System.exit(1);
        }

        if ("covariance".equals(args[0]) && args.length != 2) {
            printUsage();
            System.exit(1);
        }

        if ("correlation".equals(args[0]) && args.length != 2) {
            printUsage();
            System.exit(1);
        }

        if ("transform".equals(args[0]) && args.length != 3) {
            printUsage();
            System.exit(1);
        }

        Path path1, path2, path3;
        try {
            switch (args[0].toLowerCase()) {
                case "evd":
                    path1 = Paths.get(args[1]);
                    evd(path1);
                    break;
                case "scores":
                    path1 = Paths.get(args[1]);
                    path2 = Paths.get(args[2]);
                    path3 = null;
                    if(args.length == 4)
                    {
                    	path3 = Paths.get(args[3]);
                    }	
                    scores(path1, path2,path3);
                    break;
                case "transpose":
                    path1 = Paths.get(args[1]);
                    transpose(path1);
                    break;
                case "center":
                    path1 = Paths.get(args[1]);
                    center(path1);
                    break;
                case "scale":
                    path1 = Paths.get(args[1]);
                    scale(path1);
                    break;
                case "covariance":
                    path1 = Paths.get(args[1]);
                    covariance(path1);
                    break;
                case "correlation":
                    path1 = Paths.get(args[1]);
                    //return correlation(path1);
                case "transform":
                    path1 = Paths.get(args[1]);
                    path2 = Paths.get(args[2]);
                    transform(path1, path2);
                    break;
                default:
                    printUsage();
                    System.exit(1);
            }
        } catch (IOException ex) {
            System.err.println("Could not read or write file: " + ex.getMessage());
        } catch (NotConvergedException ex) {
            System.err.println("Eigenvector decomposition did not converge: " + ex.getMessage());
        }
		//return null;
    }

    /**
     *
     * Prints program usage instructions to stdout.
     *
     */
    private static void printUsage() {

        System.out.println("Usages:");
        System.out.println("java -jar PCA.jar center file");
        System.out.println("java -jar PCA.jar scale file");
        System.out.println("java -jar PCA.jar transpose file");
        System.out.println("java -jar PCA.jar covariance file");
        System.out.println("java -jar PCA.jar correlation file");
        System.out.println("java -jar PCA.jar evd symmetricmatrixfile");
        System.out.println("java -jar PCA.jar scores eigenvectorfile originalfile");
        System.out.println("java -jar PCA.jar transform scorefile eigenvaluefile");
    }
   
    /**
     *
     * Eigenvector decomposition.
     *
     * Calculates right eigenvectors and (real parts of) eigenvalues for a given
     * matrix.
     *
     * The given path has to contain a symmetric matrix. Uses natives and is
     * multi-threaded when natives are available. Single-threaded when natives
     * are not available.
     *
     * Writes two files: path.eigenvectors.txt and path.eigenvalues.txt
     *
     * @param path Path of matrix file to run decomposition on
     * @throws IOException If cannot read from / write to disk
     * @throws NotConvergedException If eigenvector decomposition doesn't
     * converge
     */
    public static void evd(Path path) throws IOException, NotConvergedException {
    	log("Reading data");

        String[] headers = FileUtil.readColumnHeaders(path);
        int matrixSize = headers.length;
        log("Matrix size " + matrixSize + " x " + matrixSize);

        String[] pcHeaders = new String[matrixSize];
        for (int i = 0, len = pcHeaders.length; i < len; i++) {
            pcHeaders[i] = "PC" + (i + 1);
        }

        DenseMatrix matrix = new DenseMatrix(matrixSize, matrixSize);
        FileUtil.readMatrix(path, matrix);

        log("Data read");
        MatrixStruct matrixStruct = new MatrixStruct(matrix,pcHeaders,headers);
        evd(matrixStruct, path);
    }
    
    public static MatrixStruct[] evd(MatrixStruct matrixStruct, Path path) throws IOException, NotConvergedException 
    {
    	String[] headers = matrixStruct.getColHeaders();
    	DenseMatrix matrix = matrixStruct.matrix;
    	int matrixSize = headers.length;

        log("Calculating eigenvectors");
        EVD evd = new EVD(matrixSize, false, true);
        evd.factor(matrix);
        log("Eigenvectors calculated");

        DenseMatrix rightEigenvectors = evd.getRightEigenvectors();
       
        MatrixStruct[] eigenValueDecomp = new MatrixStruct[2];
        double[] evs = evd.getRealEigenvalues();
        DenseMatrix eigen = new DenseMatrix(evs.length, 1);
        String[] pcHeaders = new String[evs.length];
        String[] colNames = new String[1];      
        for(int r = 0; r < evs.length; r++)
        {
        	eigen.set(r, 0, evs[r]);
        	pcHeaders[r] = "PC"+(r+1);
        }       
        eigenValueDecomp[0] = new MatrixStruct(rightEigenvectors, headers ,pcHeaders);
        eigenValueDecomp[0].transpose();
        
        log("Writing eigenvectors");
        eigenValueDecomp[0].write(path.toString().replace(".gz", "").replace(".txt", "") + ".eigenvectors.txt");
        
        log("Writing eigenvalues");
        try (FileWriter fw = new FileWriter(path.toString().replace(".gz", "").replace(".txt", "") + ".eigenvalues.txt")) {
        	FileUtil.writeMatrix(fw, eigen, pcHeaders, colNames);
        }
        eigenValueDecomp[1] = new MatrixStruct(eigen, pcHeaders, colNames);
        log("Done");
        return eigenValueDecomp;
    }
    
    /**
     *
     * Principal component scores.
     *
     * Calculates principal component scores and Cronbach's alpha values based
     * on a given eigenvector matrix and original data matrix.
     *
     * Orientation of the original matrix is automatically detected.
     *
     * Writes two files: originalPath.scores.txt and originalPath.cronbachsAlpha.txt
     *
     * @param evPath Path to an eigenvector matrix where each row is an
     * eigenvector
     * @param originalPath Path to the original data matrix
     * @param averagesPath Path to averages to be used for centering
     * If null calculates averages from <originalPath> file
     * @return 
     * @throws IOException If cannot read from / write to disk
     */
    public static void scores(Path evPath, Path originalPath) throws IOException {
    	scores(evPath, originalPath, null);
    }
    public static MatrixStruct[] scores(MatrixStruct eigenVectors, MatrixStruct expression, String originalPath) throws IOException 
    {
    	return scores(eigenVectors, expression, originalPath, false, true);
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
	private static void correctForRowAverages(MatrixStruct expression, int numRowsOri, int numColsOri, DenseMatrix originalMatrix, String originalPath, String[] rowHeaders) throws IOException 
	{
    	log("Centering each row of data");
    	//read averages if supplied
        //correct data for average
        for (int r = 0; r < numRowsOri; r++) {
            double mean = getMean(originalMatrix, r);
            for (int c = 0; c < numColsOri; c++) {
            	originalMatrix.add(r, c, -mean);//do the actual correction
            }
        }
        //write averages
        String[] colHeadersAvg = new String[]{"Averages"};
        log("Rows centered");        
	}

	private static void putInCorrectOrientation(MatrixStruct eigenVectors, MatrixStruct expression, boolean rotateBack) 
	{
		if(rotateBack)
    	{
    		if(!eigenVectors.getColHeaders()[0].equals("PC1") && eigenVectors.getRowHeaders()[0].equals("PC1") )
	    	{
	    		log("Principal components appear to be on rows and should be on colums when rotating back; transposing");
	    		eigenVectors.transpose();
	    	}
    	}
    	else
    	{
	    	if(eigenVectors.getColHeaders()[0].equals("PC1") && !eigenVectors.getRowHeaders()[0].equals("PC1") )
	    	{
	    		log("Principal components appear to be on columns and should be on rows; transposing");
	    		eigenVectors.transpose();
	    	}
    	}
    	if(expression.getRowHeaders().length != eigenVectors.getColHeaders().length)
    	{
    		expression.transpose();
    		log("Dimensions are incorrect, transposing expression matrix");
    		if(expression.getRowHeaders().length != eigenVectors.getColHeaders().length)
    		{
    			System.out.println("expression.rowHeaders " + expression.getRowHeaders().length + " eigenVectors.colHeaders " + eigenVectors.getColHeaders().length);
    			log("Dimensions are still incorrect, exiting");
    			System.exit(0);
    		}
    	}
		//System.out.println(x);
	}

	public static void scores(Path evPath, Path originalPath, Path averagesPath) throws IOException {

        log("Reading eigenvectors");
        String[] EVrowNames = FileUtil.readRowHeaders(evPath);
        String[] EVcolNames = FileUtil.readColumnHeaders(evPath);
        int numRowsEV = EVrowNames.length;
        int numColsEV = EVcolNames.length;
        log("Matrix size " + numRowsEV + " x " + numColsEV);
        DenseMatrix evMatrix = new DenseMatrix(numRowsEV, numColsEV);
        FileUtil.readMatrix(evPath, evMatrix);
        log("Eigenvectors read");

        log("Reading original data");
        String[] colHeaders = FileUtil.readColumnHeaders(originalPath);
        String[] rowHeaders = FileUtil.readRowHeaders(originalPath);
        int numRowsOri = rowHeaders.length;
        int numColsOri = colHeaders.length;
        DenseMatrix originalMatrix;
        if (numColsEV != numRowsOri) {//I had to change it, otherwise it fries my brains...
            colHeaders= FileUtil.readRowHeaders(originalPath);
            rowHeaders= FileUtil.readColumnHeaders(originalPath);
            numRowsOri = rowHeaders.length;
            numColsOri = colHeaders.length;
            log("Matrix size (transposed) " + numRowsOri + " x " + numColsOri);
            originalMatrix = new DenseMatrix(numRowsOri, numColsOri);
            FileUtil.readMatrix(originalPath, originalMatrix, true);
        } else {
            log("Matrix size " + numRowsOri + " x " + numColsOri);
            originalMatrix = new DenseMatrix(numRowsOri, numColsOri);
            FileUtil.readMatrix(originalPath, originalMatrix);
        }
        log("Original data read");
        MatrixStruct expression = new MatrixStruct(originalMatrix, rowHeaders, colHeaders);
        MatrixStruct eigenVectors = new MatrixStruct(evMatrix, EVrowNames, EVcolNames);
        
        scores(eigenVectors, expression, originalPath.toString());
             
        log("Done");
    }
	/**rotates matrix back (does opposite of scores)*/
    public static MatrixStruct[] rotateBack(MatrixStruct eigenVectors, MatrixStruct expression, String originalPath) throws IOException 
    {
    	return scores(eigenVectors, expression, originalPath, true, false);
    }
    /**
     *
     * Calculates mean of a given row.
     *
     * @param matrix Data matrix
     * @param row Row index
     * @return Arithmetic mean of the row
     */
    private static double getMean(DenseMatrix matrix, int row) {
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
    private static double getVariance(DenseMatrix matrix, double mean, int row) {

        double variance = 0;
        for (int c = 0, numCols = matrix.numColumns(); c < numCols; c++) {
            variance += (matrix.get(row, c) - mean) * (matrix.get(row, c) - mean) / (numCols - 1);
        }
        return variance;
    }
    
    /**
     *
     * Center each row in a given matrix.
     *
     * For each element, subtracts the mean of the corresponding row, thus
     * centering each row to a mean of zero.
     *
     * Writes path.centered.txt
     *
     * @param path Path to the data matrix
     * @throws IOException If cannot read from / write to disk
     */
    public static void center(Path path) throws IOException {

        log("Reading data");
        String[] colHeaders = FileUtil.readColumnHeaders(path);
        String[] rowHeaders = FileUtil.readRowHeaders(path);
        int numCols = colHeaders.length;
        int numRows = rowHeaders.length;
        DenseMatrix matrix = new DenseMatrix(numRows, numCols);
        FileUtil.readMatrix(path, matrix);
        log("Data read");

        log("Centering each row of data");
        for (int r = 0; r < numRows; r++) {
            double mean = getMean(matrix, r);
            for (int c = 0; c < numCols; c++) {
                matrix.add(r, c, -mean);
            }
        }
        log("Rows centered");

        log("Writing centered data");
        try (FileWriter fw = new FileWriter(path.toString().replace(".gz", "").replace(".txt", "") + ".centered.txt")) {
            FileUtil.writeMatrix(fw, matrix, rowHeaders, colHeaders);
        }

        log("Done");
    }

    /**
     *
     * Scale each row in a given matrix.
     *
     * Divides each element with the standard deviation of the corresponding
     * row, thus scaling each row to a standard deviation of one.
     *
     * Writes path.scaled.txt
     *
     * @param path Path to the data matrix
     * @throws IOException If cannot read from / write to disk
     */
    public static void scale(Path path) throws IOException {

        log("Reading data");
        String[] colHeaders = FileUtil.readColumnHeaders(path);
        String[] rowHeaders = FileUtil.readRowHeaders(path);
        int numCols = colHeaders.length;
        int numRows = rowHeaders.length;
        DenseMatrix matrix = new DenseMatrix(numRows, numCols);
        FileUtil.readMatrix(path, matrix);
        log("Data read");

        log("Scaling each row of data to a standard deviation of one");
        for (int r = 0; r < numRows; r++) {
            double mean = getMean(matrix, r);
            double variance = getVariance(matrix, mean, r);
            double stDev = Math.sqrt(variance);
            for (int c = 0; c < numCols; c++) {
                matrix.set(r, c, matrix.get(r, c) / stDev);
            }
        }
        log("Rows scaled");

        log("Writing scaled data");
        try (FileWriter fw = new FileWriter(path.toString().replace(".gz", "").replace(".txt", "") + ".scaled.txt")) {
            FileUtil.writeMatrix(fw, matrix, rowHeaders, colHeaders);
        }

        log("Done");
    }

    /**
     *
     * Covariance over rows.
     *
     * Calculates covariance for each pair of rows in a given matrix.
     *
     * Writes the symmetric covariance matrix path.covariance.txt
     *
     * @param path Path to the data matrix
     * @throws IOException If cannot read from / write to disk
     */
    public static MatrixStruct covariance(Path path) throws IOException 
    {
    	return covariance(path, null);
    }
    public static MatrixStruct covariance(Path path, String writeName) throws IOException {

        log("Reading data");
        String[] colHeaders = FileUtil.readColumnHeaders(path);
        String[] rowHeaders = FileUtil.readRowHeaders(path);
        int numCols = colHeaders.length;
        int numRows = rowHeaders.length;
        DenseMatrix matrix = new DenseMatrix(numRows, numCols);
        FileUtil.readMatrix(path, matrix);
        log("Data read");

        DenseMatrix covarianceMatrix = new DenseMatrix(numRows, numRows);
        log("Calculating covariance matrix");
        for (int r1 = 0; r1 < numRows; r1++) {
            double mean1 = getMean(matrix, r1);
            for (int r2 = r1; r2 < numRows; r2++) {
                double mean2 = getMean(matrix, r2);
                double covariance = 0;
                for (int c = 0; c < numCols; c++) {
                    covariance += (matrix.get(r1, c) - mean1) * (matrix.get(r2, c) - mean2) / (numCols - 1);
                }
                covarianceMatrix.set(r1, r2, covariance);
                covarianceMatrix.set(r2, r1, covariance);
            }
        }
        log("Covariance matrix calculated");

        log("Writing covariance matrix");
        if(writeName == null)
        	writeName = path.toString().replace(".gz", "").replace(".txt", "") + ".covariance.txt";
        try (FileWriter fw = new FileWriter(writeName)) {
            FileUtil.writeMatrix(fw, covarianceMatrix, rowHeaders, rowHeaders);
        }

        log("Done");
        MatrixStruct covMat = new MatrixStruct(covarianceMatrix, rowHeaders, rowHeaders);
        return covMat;
    }

    /**
     *
     * Pearson correlation over rows.
     *
     * Calculates Pearson correlation for each pair of rows in a given matrix.
     *
     * Writes the symmetric correlation matrix path.correlation.txt
     *
     * @param path Path to the data matrix
     * @throws IOException If cannot read from / write to disk
     */
    public static MatrixStruct[] correlation(Path path) throws IOException {

        log("Reading data");
        String[] colHeaders = FileUtil.readColumnHeaders(path);
        String[] rowHeaders = FileUtil.readRowHeaders(path);
        int numCols = colHeaders.length;
        int numRows = rowHeaders.length;
        DenseMatrix matrix = new DenseMatrix(numRows, numCols);
        FileUtil.readMatrix(path, matrix);
        log("Data read");

        DenseMatrix correlationMatrix = new DenseMatrix(numRows, numRows);
        log("Calculating correlation matrix");
        for (int r1 = 0; r1 < numRows; r1++) {
            correlationMatrix.set(r1, r1, 1);
            double mean1 = getMean(matrix, r1);
            double var1 = getVariance(matrix, mean1, r1);
            for (int r2 = r1 + 1; r2 < numRows; r2++) {
                double mean2 = getMean(matrix, r2);
                double var2 = getVariance(matrix, mean2, r2);
                double covariance = 0;
                for (int c = 0; c < numCols; c++) {
                    covariance += (matrix.get(r1, c) - mean1) * (matrix.get(r2, c) - mean2) / (numCols - 1);
                }
                double denom = Math.sqrt(var1 * var2);
                correlationMatrix.set(r1, r2, covariance / denom);
                correlationMatrix.set(r2, r1, covariance / denom);
            }
        }
        log("Correlation matrix calculated");

        log("Writing correlation matrix");
        String writeName =path.toString().replace(".gz", "").replace(".txt", "") + ".correlation.txt";
        try (FileWriter fw = new FileWriter(writeName)) {
            FileUtil.writeMatrix(fw, correlationMatrix, rowHeaders, rowHeaders);
        }

        log("Done");
        
        MatrixStruct[] output = new MatrixStruct[1];
        output[0] = new MatrixStruct(correlationMatrix, rowHeaders, rowHeaders);
        return output;
    }

    /**
     *
     * Transpose a given data matrix.
     *
     * Writes path.transposed.txt
     *
     * @param path Path to the data matrix
     * @throws IOException If cannot read from / write to disk
     */
    public static void transpose(Path path) throws IOException {

        log("Reading data");
        String[] colHeaders = FileUtil.readColumnHeaders(path);
        String[] rowHeaders = FileUtil.readRowHeaders(path);
        int numCols = colHeaders.length;
        int numRows = rowHeaders.length;
        DenseMatrix matrix = new DenseMatrix(numRows, numCols);
        FileUtil.readMatrix(path, matrix);
        log("Data read");

        log("Writing transposed data");
        try (FileWriter fw = new FileWriter(path.toString().replace(".gz", "").replace(".txt", "") + ".transposed.txt")) {
            FileUtil.writeMatrix(fw, matrix, colHeaders, rowHeaders, true);
        }

        log("Done");
    }

    /**
     *
     * Change orientation of PCA results.
     *
     * Calculates "transposed eigenvectors" from a given principal component score
     * file and eigenvalue file.
     *
     * Note: This works in the following scenario:
     *
     * - There is an original input matrix X
     *
     * - Each column of X has been centered (zero mean for each column)
     *
     * - A covariance matrix has been calculated over rows (for pairs of rows) of X
     * 
     * - Eigenvector decomposition has been done for this covariance matrix
     * 
     * - Principal component scores have been calculated based on the resulting eigenvectors
     *
     * Then, this method calculates "transposed eigenvectors": eigenvectors as
     * they would be if PCA had been done this way:
     * 
     * - Original input matrix Y = X' (X transposed)
     * 
     * - Each column of Y has been centered (zero mean for each column)
     * 
     * - A covariance matrix has been calculated over rows (for pairs of rows) of Y
     * 
     * - Eigenvector decomposition has been done for this covariance matrix
     * 
     * Writes scorePath.transformedToEigenvectors.txt
     * 
     * @param scorePath Path to a principal component score matrix
     * where each row is a component
     * @param eigenvaluePath Path to a file with eigenvalues
     * (corresponding to the eigenvectors from which the principal component scores have been calculated)
     * @return 
     * @throws IOException If cannot read from / write to disk
     */
    public static MatrixStruct transform(MatrixStruct scores, MatrixStruct eigenValues, String scorePath) throws IOException 
    {
    	
    	String[] rowHeaders = scores.getRowHeaders();
        String[] colHeaders = scores.getColHeaders();
        int numCols = colHeaders.length;
        int numRows = rowHeaders.length;
        DenseMatrix scoreMatrix = scores.matrix;
        DenseMatrix eigenvalueMatrix = eigenValues.matrix;
        
        log("Calculating transposed eigenvectors");
        DenseMatrix transposedEVMatrix = new DenseMatrix(numRows, numCols);
        for (int comp = 0; comp < numRows; comp++) {
            for (int i = 0; i < numCols; i++) {
                double score = scoreMatrix.get(comp, i);
                Double newScore = score / (Math.sqrt((numCols-1) * eigenvalueMatrix.get(comp, 0)));//I(Sipko) added -1 to numCols
                if(newScore.isNaN())
                	newScore = (double) 0;
                transposedEVMatrix.set(comp, i, newScore);
            }
        }
        log("Transformation calculated");
        String writeFN = scorePath;
        try (FileWriter fw = new FileWriter(writeFN)) {
            FileUtil.writeMatrix(fw, transposedEVMatrix, rowHeaders, colHeaders);
        }
        MatrixStruct transposed = new MatrixStruct(transposedEVMatrix, rowHeaders, colHeaders);
        log("Done");
        return transposed;
    }
    public static void transform(Path scorePath, Path eigenvectorPath) throws IOException {

        log("Reading scores");
        String[] colHeaders = FileUtil.readColumnHeaders(scorePath);
        String[] rowHeaders = FileUtil.readRowHeaders(scorePath);
        int numCols = colHeaders.length;
        int numRows = rowHeaders.length;
        DenseMatrix scoreMatrix = new DenseMatrix(numRows, numCols);
        FileUtil.readMatrix(scorePath, scoreMatrix);
        log("Scores read");

        log("Reading eigenvalues");
        DenseMatrix eigenvalueMatrix = new DenseMatrix(numRows, 1);
        FileUtil.readArray(eigenvectorPath, eigenvalueMatrix, 1);
        log("Eigenvalues read");
        
        MatrixStruct scores = new MatrixStruct(scoreMatrix, colHeaders,rowHeaders);
        MatrixStruct eigenvalue = new MatrixStruct(eigenvalueMatrix); 
        transform(scores, eigenvalue, scorePath.toString().replace(".gz", "").replace(".txt", "") + ".transformedToEigenvectors.txt");
    }

    /**
     * 
     * Calculate Cronbach's alpha for each principal component.
     * 
     * @param evMatrix Eigenvector matrix, each row is an eigenvector
     * @param scoreMatrix Principal component score matrix, each row is a component
     * @return An array of Cronbach's alpha values
     */
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

    /**
     * 
     * Print given text to stdout, preceded by a MySQL-style timestamp and a tab.
     * 
     * @param text Text to log
     */
    public static void log(String text) {

        String time = timeFormat.format(new Date());
        System.out.println(time + "\t" + text);
    }
}
