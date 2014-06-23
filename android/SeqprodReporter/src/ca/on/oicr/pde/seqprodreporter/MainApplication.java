package ca.on.oicr.pde.seqprodreporter;

import android.app.Application;

public class MainApplication extends Application {
	private boolean isCurrentActivityVisible;
	
	public boolean getisCurrentActivityVisible(){
		return isCurrentActivityVisible;
	}

	public void setisCurrentActivityVisible(boolean isVisible){
		isCurrentActivityVisible = isVisible;
	}
}
