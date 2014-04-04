package ca.on.oicr.pde.seqprodreporter;

import android.os.Bundle;
import android.preference.Preference;
import android.preference.PreferenceFragment;
import android.preference.PreferenceScreen;

public class SeqprodPreferencesFragment extends PreferenceFragment 
    implements Preference.OnPreferenceChangeListener {
		
	public SeqprodPreferencesFragment() {
		super();
	}
	
	public void onCreate(Bundle savedInstanceState) {
	super.onCreate(savedInstanceState);

	try {
	    addPreferencesFromResource(R.xml.preferences);
	    PreferenceScreen ps = getPreferenceScreen();
	    for (int i = 0; i < ps.getPreferenceCount(); i++) {
	      ps.getPreference(i).setOnPreferenceChangeListener(this);
	    }
	} catch (Exception e) {
		//TODO there maybe a conflict with stored SharedPreferences values, 
		//in this case delete the stored values and try again
		e.printStackTrace();
	}

    }

	@Override
	public boolean onPreferenceChange(Preference preference, Object newValue) {
		// TODO Auto-generated method stub
		Preference p = getPreferenceScreen().findPreference(preference.getKey());
		if (null != p)
		 p.setSummary("Changed");
		return true;
	}

}
