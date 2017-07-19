package PizzlyClasses;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import Tools.Script;

import Tools.Script;

public class PizzlyFusionStructure extends Script<PizzlyFusionStructure>
{
	Fusion[] genes = null;

	public Fusion[] getFusions()
	{
		return genes;
	}

	public void getFusions(Fusion[] genes)
	{
		this.genes = genes;
	} 
}
