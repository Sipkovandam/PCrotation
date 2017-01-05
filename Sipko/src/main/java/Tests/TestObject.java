package Tests;

import java.util.HashMap;

import Tools.Script;

public class TestObject extends Script implements Cloneable
{

	int f =1;
	String bla = "Bla";
	String[] ts = null;
	HashMap<String,String> hash= null;
	
	public Object clone() throws CloneNotSupportedException {
	      return super.clone();
	  }
	
	public void run()
	{
		p("testOut");
	}
}
