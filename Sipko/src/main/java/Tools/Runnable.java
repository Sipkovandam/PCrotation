package Tools;

import java.util.HashMap;
import java.util.List;

import MatrixScripts.MyMatrix;

public interface Runnable {
	
	Exception notRunnableException = new Exception("Function \"run()\" not declared");
	Exception notRunnableExceptionWithArg = new Exception("Function \"run(Object var)\" not declared");
	Exception notRunnableExceptionWithArgs = new Exception("Function \"run(Object ... var)\" not declared");
	Exception notRunnableExceptionWithMyMatrix= new Exception("Function \"run(MyMatrix passMatrix)\" not declared");

	default void run() throws Exception
	{
		throw notRunnableException;
	}

	
	default Object run(Object ... var) throws Exception
	{
		throw notRunnableExceptionWithArgs;
	}
	
	default HashMap<String,Object> run(HashMap<String,Object> vars) throws Exception
	{
		throw notRunnableException;
	}


	default Object run(Object var) throws Exception
	{
		throw notRunnableExceptionWithArg;
	}
	
	default MyMatrix run(MyMatrix passMatrix) throws Exception
	{
		this.run();
		return null;
	}
}
