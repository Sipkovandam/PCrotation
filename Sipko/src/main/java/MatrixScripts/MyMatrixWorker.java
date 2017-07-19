package MatrixScripts;

public class MyMatrixWorker implements Runnable
{
	private MyMatrix matrix = null;
	private String fileName = null;
	
	public MyMatrixWorker(MyMatrix matrix, String fileName)
	{
		this.matrix = matrix;
		this.fileName = fileName;
	}
	
	@Override
	public void run()
	{
		matrix.write(fileName);
	}
}
