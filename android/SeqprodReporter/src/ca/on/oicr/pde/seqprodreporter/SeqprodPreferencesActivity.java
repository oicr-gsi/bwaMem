package ca.on.oicr.pde.seqprodreporter;

import android.app.Activity;
import android.os.Bundle;

/**
 * The Activity that contains the SeqprodPreferencesFragment
 *
 * @see SeqprodPreferencesFragment
 */
public class SeqprodPreferencesActivity extends Activity {
	
	public void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		//startWithFragment("Preferences",savedInstanceState, new SeqprodPreferencesFragment(), 0);
		setContentView(R.layout.pref_fragment_container);
		
		getFragmentManager().beginTransaction()
		 .replace(android.R.id.content, new SeqprodPreferencesFragment())
		 .commit();
	}
}
