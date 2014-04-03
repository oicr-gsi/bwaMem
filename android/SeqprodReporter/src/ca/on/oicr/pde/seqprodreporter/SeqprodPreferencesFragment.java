package ca.on.oicr.pde.seqprodreporter;

import android.os.Bundle;
import android.preference.PreferenceFragment;

public class SeqprodPreferencesFragment extends PreferenceFragment {
		
	public SeqprodPreferencesFragment() {
		super();
	}
	
	public void onCreate(Bundle savedInstanceState) {
	super.onCreate(savedInstanceState);

	try {
	    addPreferencesFromResource(R.xml.preferences);
	} catch (Exception e) {
		//TODO there maybe a conflict with stored SharedPreferences values, 
		//in this case delete the stored values and try again
		e.printStackTrace();
	}

    }

}
