package Tests;

import Tools.Script;

public class TestObject3 <M> extends Script<TestObject3<M>> 
{
	String iets = "iets";
	TestObject testObject = new TestObject("bla");
	TestObject2 testObject2 = new TestObject2("blabla", 15.0);
	M testObject3 = null;

	public TestObject3()
	{
		testObject3 = (M) new TestObject("zzzz");
		TestObject test = (TestObject) testObject3;
		test.log(test.bla);
		setJsonFN("E:/Groningen/Test/JSON/TestCombinedObjects.config");
	}
	
	public void run()
	{
		testObject3 = null;
		log("ruunning!");
	}
}
