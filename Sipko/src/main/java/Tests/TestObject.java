package Tests;

import java.util.HashMap;

import Tools.Script;

public class TestObject extends Script<TestObject> implements Cloneable, Runnable
{

	int f =1;
	String bla = "Bla";
	String[] ts = new String[]{"akkel","kaas","pizza","friet"};
	//HashMap<String,String> hash= null;
	public TestObject(String bla)
	{
		this.bla = bla;
	}
	
	public Object clone() throws CloneNotSupportedException 
	{
	      return super.clone();
	}
	
	public void run()
	{
		p("Running!");
	}
}
