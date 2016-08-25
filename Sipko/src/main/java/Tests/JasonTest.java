package Tests;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonObject;

import PCA.Var;

public class JasonTest 
{
	
	static Var var = null;
	public static void main(String[] args)
	{
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
	      
	      var = new Var();
	      var.jsonFN = "E:/Groningen/Data/Test/test.txt";
	      var.writeVars();
	      var.readVars("E:/Groningen/Data/Test/test.txt");
	 }
}
