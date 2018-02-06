package Listeners.FileNameListeners;

import Listeners.Listener;

public class FileNameListener extends Listener
{
	@Override
	public void run()
	{
		log("FileNameListener requires fileName as input variable");
	}
	
	public void run(String fileName)
	{
		log("Run(String) not defined for Listener");
	}
}
