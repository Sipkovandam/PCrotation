package Listeners.FileNameListeners;

import Listeners.Listener;

public class FileNameListener extends Listener
{
	@Override
	public void run()
	{
		p("FileNameListener requires fileName as input variable");
	}
	
	public void run(String fileName)
	{
		p("Run(String) not defined for Listener");
	}
}
