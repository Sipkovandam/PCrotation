package PCA;

import java.io.IOException;

import Tools.FileUtils;
import Tools.JSONutil;

public class RlogLargeMatrix_Main 
{

	public static void main(String[] args) throws IOException
	{
		//args = new String[]{"json=E:/Groningen/Test/JSON/jsonTest.json"};
		RlogLargeMatrix rlogLargeMatrix = new JSONutil<RlogLargeMatrix>().read(args, new RlogLargeMatrix());
		rlogLargeMatrix=rlogLargeMatrix.checkArgs(args);
		rlogLargeMatrix.run();
	}
}
