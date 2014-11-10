package ca.on.oicr.pde.seqprodreporter;

import android.app.AlertDialog;
import android.content.DialogInterface;
import android.content.DialogInterface.OnClickListener;
import android.content.Intent;
import android.content.SharedPreferences;
import android.os.Bundle;
import android.preference.Preference;
import android.preference.PreferenceFragment;
import android.support.v4.content.LocalBroadcastManager;
import android.webkit.URLUtil;

/**
 * An extended PreferenceFragment that is used to display the preferences that a user will set.
 *
 * @see PreferenceFragment
 *
 */
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
		getPreferenceManager().getSharedPreferences().unregisterOnSharedPreferenceChangeListener(this);
		super.onStop();
	}

	public SeqprodPreferencesFragment() {
		super();
	}
	
	public void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		getPreferenceManager().setSharedPreferencesName(
				ReporterActivity.PREFERENCE_FILE);

		try {
			addPreferencesFromResource(R.xml.preferences);
			SharedPreferences sp = getPreferenceManager()
					.getSharedPreferences();
			for (int s = 0; s < getPreferenceScreen().getPreferenceCount(); s++) {
				Preference p = getPreferenceScreen().getPreference(s);
				if (null != p) {
					String existingValue = sp.getString(p.getKey(), "NA");
					if (!existingValue.equals("NA"))
						p.setSummary(existingValue);
				}
			}
		} catch (Exception e) {
			// TODO there maybe a conflict with stored SharedPreferences values,
			// in this case delete the stored values and try again
			e.printStackTrace();
		}

	}

	/*
	 * (non-Javadoc) * @see
	 * android.content.SharedPreferences.OnSharedPreferenceChangeListener
	 * #onSharedPreferenceChanged(android.content.SharedPreferences,
	 * java.lang.String) Here we update summary text inside preference fragment
	 * and send a broadcast (to the receiver in ReporterActivity)
	 */
	@Override
	public void onSharedPreferenceChanged(SharedPreferences sharedPreferences,
			String key) {
		if (key.equals("pref_hostName")) {
			String newUrl = getPreferenceManager().getSharedPreferences().getString(key, "NA");
			if (null != newUrl && !newUrl.isEmpty() && !newUrl.equals("NA")) {
				if (!URLUtil.isValidUrl(newUrl)) {

					AlertDialog.Builder builder = new AlertDialog.Builder(
							getActivity());
					builder.setMessage(R.string.url_invalid).setTitle("Error");
					builder.setNegativeButton("Close", new OnClickListener() {

						@Override
						public void onClick(DialogInterface dialog, int which) {
							getPreferenceManager().getSharedPreferences()
									.edit().remove("pref_hostName").commit();
						}
					});
					builder.create().show();
					getPreferenceScreen().findPreference("pref_hostName")
							.setSummary(R.string.pref_hostName_default);
					return;
				}
			} else { // Disable other buttons
				getPreferenceScreen().findPreference("pref_summaryScope")
						.setEnabled(false);
				getPreferenceScreen().findPreference("pref_syncFreq")
						.setEnabled(false);
			}
		}

		this.prefsChanged = true;
		Preference p = getPreferenceScreen().findPreference(key);

		if (null != p) {
			String updatedValue = getPreferenceManager().getSharedPreferences()
					.getString(key, "NA");
			if (null != updatedValue && !updatedValue.equals("NA")) {
				p.setSummary(updatedValue);
				if (key.equals("pref_hostName")) {
					getPreferenceScreen().findPreference("pref_summaryScope")
							.setEnabled(true);
					getPreferenceScreen().findPreference("pref_syncFreq")
							.setEnabled(true);
				}
			}

		}
	}
}
