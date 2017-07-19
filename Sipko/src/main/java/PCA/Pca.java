package PCA;

import java.io.IOException;
import java.nio.file.Paths;

import MatrixScripts.MyMatrix;
import Tools.Script;
import no.uib.cipr.matrix.NotConvergedException;

public class Pca extends Script<Pca>
{
	//conducts PCA over input matrix using Java. This function is slow and it is advised to use it only for matrixes with less then 10,000 rows
	String inputMatrixComment = "/root/directory/matrix.txt; MANDATORY; Matrix over which the PCA should be ran";
	String inputMatrix = null;
	String eigenVectorWriteFnComment = "/root/directory/matrix.txt; OPTIONAL; Filename to which the eigenvector file should be written";
	String eigenVectorWriteFn = null;
	String eigenValueWriteComment = "/root/directory/matrix.txt; OPTIONAL;Filename to which the eigenvalue file should be written";
	String eigenValueWriteFn = null;
	
	@Override
	public void run()
	{
		try
		{
			JuhaPCA.PCA.evd(Paths.get(inputMatrix), eigenVectorWriteFn, eigenValueWriteFn);
		} catch (Exception e)
		{
			e.printStackTrace();
		}
	}
	
	@Override
	public MyMatrix run(MyMatrix input)
	{
		try
		{
			JuhaPCA.PCA.evd(Paths.get(inputMatrix), eigenVectorWriteFn, eigenValueWriteFn);
		} catch (Exception e)
		{
			e.printStackTrace();
		}
		return input;
	}

	public String getInputMatrix()
	{
		return inputMatrix;
	}

	public void setInputMatrix(String inputMatrix)
	{
		this.inputMatrix = inputMatrix;
	}

	public String getEigenVectorWriteFn()
	{
		return eigenVectorWriteFn;
	}

	public void setEigenVectorWriteFn(String eigenVectorWriteFn)
	{
		this.eigenVectorWriteFn = eigenVectorWriteFn;
	}

	public String getEigenValueWriteFn()
	{
		return eigenValueWriteFn;
	}

	public void setEigenValueWriteFn(String eigenValueWriteFn)
	{
		this.eigenValueWriteFn = eigenValueWriteFn;
	}
}
