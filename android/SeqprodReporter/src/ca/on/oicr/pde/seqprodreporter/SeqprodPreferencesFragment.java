package ca.on.oicr.pde.seqprodreporter;

import android.content.Intent;
import android.content.SharedPreferences;
import android.os.Bundle;
import android.preference.Preference;
import android.preference.PreferenceFragment;
import android.support.v4.content.LocalBroadcastManager;

public class SeqprodPreferencesFragment extends PreferenceFragment 
    implements SharedPreferences.OnSharedPreferenceChangeListener {
	private boolean prefsChanged;
	
	@Override
	public void onStart() {
		super.onStart();
		this.prefsChanged = false;
		getPreferenceManager().getSharedPreferences().registerOnSharedPreferenceChangeListener(this);
	}

	@Override
	public void onStop() {
		if (this.prefsChanged)
			LocalBroadcastManager.getInstance(getActivity().getApplicationContext())
			 .sendBroadcast(new Intent(ReporterActivity.PREFCHANGE_INTENT));
			//getActivity().getApplicationContext().sendBroadcast(new Intent(ReporterActivity.PREFCHANGE_INTENT));
		getPreferenceManager().getSharedPreferences().unregisterOnSharedPreferenceChangeListener(this);
		super.onStop();
	}

	public SeqprodPreferencesFragment() {
		super();
	}
	
	public void onCreate(Bundle savedInstanceState) {
	super.onCreate(savedInstanceState);
	getPreferenceManager().setSharedPreferencesName(ReporterActivity.PREFERENCE_FILE);
	
	try {
	    addPreferencesFromResource(R.xml.preferences); 
	    SharedPreferences sp = getPreferenceManager().getSharedPreferences();
	    for (int s = 0; s < getPreferenceScreen().getPreferenceCount(); s++) {
	       Preference p = getPreferenceScreen().getPreference(s);
	       if (null != p) {
	         String existingValue = sp.getString(p.getKey(), "NA");
	         if (!existingValue.equals("NA"))
	            p.setSummary(existingValue);
	       }	
	    }
	} catch (Exception e) {
		//TODO there maybe a conflict with stored SharedPreferences values, 
		//in this case delete the stored values and try again
		e.printStackTrace();
	}

    }

	/*
	 * (non-Javadoc)
	 * * @see android.content.SharedPreferences.OnSharedPreferenceChangeListener#onSharedPreferenceChanged(android.content.SharedPreferences, java.lang.String)
	 * Here we update summay text inside preference fragment and send a broadcast (to the receiver in ReporterActivity)
	 */
	@Override
	public void onSharedPreferenceChanged(SharedPreferences sharedPreferences,
			String key) {
		this.prefsChanged = true;
		Preference p = getPreferenceScreen().findPreference(key);
		if (null != p) {
		 String updatedValue = getPreferenceManager().getSharedPreferences().getString(key, "NA");
		 if (!updatedValue.equals("NA"))
			 p.setSummary(updatedValue);
		}
	}
}
