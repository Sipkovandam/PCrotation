package Tests;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonObject;

import PCA.RlogLargeMatrix;
import PCA.Var;
import Tools.JSONutil;
import Tools.Script;

public class JsonTest extends Script<JsonTest> 
{
	
	static JsonVar var = null;
	public static void main(String[] args)
	{
		var = new JsonVar();
		var.setJsonFN("E:/Groningen/Test/JSON/test.config");
		var.writeConfig();
//	      JsonObject obj = new JsonObject();
//	      GsonBuilder gsonBuilder = new GsonBuilder();
//	      Gson gson = gsonBuilder.create();
//	      
//	      obj.addProperty("name", "foo");
//	      obj.addProperty("num", new Integer(100));
//	      obj.addProperty("balance", new Double(1000.21));
//	      obj.addProperty("is_vip", new Boolean(true));
//	
//	      System.out.print(obj);
//	      System.out.print(obj.get("num").getAsInt());
//	      gson.fromJson(obj, null);
		  JsonVar rlogLargeMatrix = new JSONutil<JsonVar>().read("E:/Groningen/Data/Test/testJson.txt", new JsonVar());
//	      var = new JsonVar();
//	      var.write("E:/Groningen/Data/Test/testJson.txt", var);
	      System.out.println("done");
//	      var.jsonFN = "E:/Groningen/Data/Test/test.txt";
//	      var.writeVars();
//	      var.readVars("E:/Groningen/Data/Test/test.txt");
	 }
}
