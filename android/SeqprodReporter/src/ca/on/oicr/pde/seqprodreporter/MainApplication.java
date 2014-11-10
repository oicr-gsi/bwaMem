package ca.on.oicr.pde.seqprodreporter;

import android.app.Application;

/**
 * Used by the application's notification system to see whether or not the application
 * is in the foreground or not.
 * 
 */
public class MainApplication extends Application {
	private boolean isCurrentActivityVisible;
	
	public boolean getisCurrentActivityVisible(){
		return isCurrentActivityVisible;
	}

	public void setisCurrentActivityVisible(boolean isVisible){
		isCurrentActivityVisible = isVisible;
	}
}
