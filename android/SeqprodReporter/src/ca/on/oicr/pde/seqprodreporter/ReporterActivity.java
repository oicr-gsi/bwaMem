package ca.on.oicr.pde.seqprodreporter;

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Timer;
import java.util.TimerTask;

import android.content.BroadcastReceiver;
import android.content.Context;
import android.content.Intent;
import android.content.IntentFilter;
import android.content.SharedPreferences;
import android.os.Bundle;
import android.support.v4.app.Fragment;
import android.support.v4.app.FragmentManager;
import android.support.v4.app.FragmentPagerAdapter;
import android.support.v4.app.FragmentTransaction;
import android.support.v4.content.LocalBroadcastManager;
import android.support.v4.view.ViewPager;
import android.support.v7.app.ActionBar;
import android.support.v7.app.ActionBarActivity;
import android.util.Log;
import android.view.Menu;
import android.view.MenuItem;

public class ReporterActivity extends ActionBarActivity implements
		ActionBar.TabListener {

	/**
	 * The {@link android.support.v4.view.PagerAdapter} that will provide
	 * fragments for each of the sections. We use a {@link FragmentPagerAdapter}
	 * derivative, which will keep every loaded fragment in memory. If this
	 * becomes too memory intensive, it may be best to switch to a
	 * {@link android.support.v4.app.FragmentStatePagerAdapter}.
	 */
	SectionsPagerAdapter mSectionsPagerAdapter;
	private final PreferenceUpdateReceiver prefUpdateReceiver = new PreferenceUpdateReceiver();

	/**
	 * The {@link ViewPager} that will host the section contents.
	 */
	ViewPager mViewPager;
    
    protected static String[] types = {"completed","failed","pending"};
    protected static String SYNC_OFF;
    public static final String TAG = "Seqprodbio Reporter";
    public static final String TIMER_NAME = "Seqprodbio Timer";
    public static final String PREFERENCE_FILE = "seqprod.conf";
    public static final String DATA_FILE = "seqprod_data.json";
    private static final long INITIAL_TIMER_DELAY = 5 * 1000L;
    
    static final String PREFCHANGE_INTENT = "ca.on.oicr.pde.seqprodreporter.prefsChanged";
      
    private String updateHost;
    private String updateRange;
    private int updateFrequency; // in minutes
    private boolean timerScheduled = false;
    
    private SharedPreferences sp;
    private Timer timer;

	@Override
	protected void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		SYNC_OFF = getResources().getString(
				R.string.pref_automaticUpdates_default);
		
		setContentView(R.layout.activity_reporter);
		//Set up timer for regular updates
		//this.timer = new Timer(TIMER_NAME);
		//Register receiver for preference updates
		LocalBroadcastManager lmb = LocalBroadcastManager.getInstance(this);
		IntentFilter prefchangeFilter = new IntentFilter(PREFCHANGE_INTENT);
		lmb.registerReceiver(prefUpdateReceiver, prefchangeFilter);
		
		this.sp = getSharedPreferences(PREFERENCE_FILE, MODE_PRIVATE);
		//Read preferences
		this.updateActivityPrefs();
		
		// Set up the action bar.
		final ActionBar actionBar = getSupportActionBar();
		actionBar.setNavigationMode(ActionBar.NAVIGATION_MODE_TABS);
        
		// Create the adapter that will return a fragment for each of the three
		// primary sections of the activity.
		mSectionsPagerAdapter = new SectionsPagerAdapter(
				getSupportFragmentManager());

		// Set up the ViewPager with the sections adapter.
		mViewPager = (ViewPager) findViewById(R.id.pager);
		mViewPager.setAdapter(mSectionsPagerAdapter);

		// When swiping between different sections, select the corresponding
		// tab. We can also use ActionBar.Tab#select() to do this if we have
		// a reference to the Tab.
		mViewPager
				.setOnPageChangeListener(new ViewPager.SimpleOnPageChangeListener() {
					@Override
					public void onPageSelected(int position) {
						actionBar.setSelectedNavigationItem(position);
					}
				});

		// For each of the sections in the app, add a tab to the action bar.
		for (int i = 0; i < mSectionsPagerAdapter.getCount(); i++) {
			// Create a tab with text corresponding to the page title defined by
			// the adapter. Also specify this Activity object, which implements
			// the TabListener interface, as the callback (listener) for when
			// this tab is selected.
			actionBar.addTab(actionBar.newTab()
					.setText(mSectionsPagerAdapter.getPageTitle(i))
					.setTabListener(this));
		}
	}
	
	@Override
	protected void onStop() {
		// TODO Switch From showing toaster messages to Notifications
		super.onStop();
	}
	
	@Override
	public void onDestroy() {
		this.timer.cancel();
		super.onDestroy();
	}
	

	@Override
	public boolean onCreateOptionsMenu(Menu menu) {
		// Inflate the menu; this adds items to the action bar if it is present.
		getMenuInflater().inflate(R.menu.reporter, menu);
		return true;
	}

	
	@Override
	public boolean onOptionsItemSelected(MenuItem item) {
		// Handle action bar item clicks here. The action bar will
		// automatically handle clicks on the Home/Up button, so long
		// as you specify a parent activity in AndroidManifest.xml.
		int id = item.getItemId();
		if (id == R.id.action_settings) {
			Intent setPrefs = new Intent(this, SeqprodPreferencesActivity.class);
			startActivity(setPrefs);
			return true;			
		}
		return super.onOptionsItemSelected(item);
	}

	@Override
	public void onTabSelected(ActionBar.Tab tab,
			FragmentTransaction fragmentTransaction) {
		// When the given tab is selected, switch to the corresponding page in
		// the ViewPager.
		mViewPager.setCurrentItem(tab.getPosition());
	}

	@Override
	public void onTabUnselected(ActionBar.Tab tab,
			FragmentTransaction fragmentTransaction) {
	}

	@Override
	public void onTabReselected(ActionBar.Tab tab,
			FragmentTransaction fragmentTransaction) {
	}

	/**
	 * A {@link FragmentPagerAdapter} that returns a fragment corresponding to
	 * one of the sections/tabs/pages.
	 */
	public class SectionsPagerAdapter extends FragmentPagerAdapter {
		private List<ReportListFragment> fragments;
		
		public SectionsPagerAdapter(FragmentManager fm) {
			super(fm);
			fragments = new ArrayList<ReportListFragment>();
		}

		@Override
		public Fragment getItem(int position) {
			// getItem is called to instantiate the fragment for the given page.
			if (position >= this.fragments.size() || null == this.fragments.get(position)) { 
			  fragments.add(position, ReportListFragment.newInstance(position + 1));
			}
			return fragments.get(position);
		}

		@Override
		public int getCount() {
			return types.length;
		}

		@Override
		public CharSequence getPageTitle(int position) {
			Locale l = Locale.getDefault();
			if (position >= 0 && position < types.length)
				return types[position].toUpperCase(l);
			return null;
		}
	}

	
	/*
	 * This function schedules the Alarm, but DOES NOT launches Http task
	 * Broadcast Receiver for Alarm will launch the actual update
	 */
	private void scheduleUpdate() {
		//TODO read shared preferences, May check validity of server Url here
		// Update data range here as well
		Log.v(TAG, "Entered sheduleUpdates");
		//Start update ONLY if host URL is valid //TODO remove Url checking
		if (null == this.updateRange || null == this.updateHost 
			|| this.updateRange.equals(getResources().getString(R.string.pref_summaryScope_default))
			|| this.updateFrequency == 0) {
			if (null != this.timer) {
			    this.timer.cancel();
			    Log.d(TAG, "Http Updates Canceled");
			}
			this.timerScheduled = false;
			this.updateFrequency = 0;
			return;
		}
		Log.d(TAG, "Will schedule Timer to trigger updates from " + this.updateHost);
		//Schedule Alert here
		//DEBUG ONLY
		this.updateFrequency = 1;
		long INTERVAL = this.updateFrequency * 15 * 1000L;
		if (this.timerScheduled) {
			this.timer.cancel();
		    Log.d(TAG, "Http Updates Canceled");	
		}
		this.timer = new Timer();
		this.timer.schedule(new TimedHttpTask(), INITIAL_TIMER_DELAY, INTERVAL);
		this.timerScheduled = true;
	}
	
	/*
	 * Update Activity's variables with new values from SharedPreferences,
	 * launch Http task if going from OFF to ON state
	 */
	private void updateActivityPrefs() {
		Log.v(TAG, "Entered updateActivityPrefs");
		if (null == this.sp) // This check may be not needed
			return;
		this.updateHost = sp.getString("pref_hostName", null);
		this.updateRange = sp.getString("pref_summaryScope", null);
		
		String uf = sp.getString("pref_syncFreq", SYNC_OFF);
		if (!uf.equals(SYNC_OFF)) {
			try {
				this.updateFrequency = Integer.parseInt(uf.substring(0, uf.lastIndexOf(" ")));
			} catch (NumberFormatException ne) {
				Log.e(TAG, "Malformed value for update frequency, will not update");
				this.updateFrequency = 0;
			}
		} else {
			this.updateFrequency = 0;
		}
		scheduleUpdate();
	}
	
	/*
	 * Broadcast Receiver for Preference Update Broadcast, updates variables with preference values
	 */
	class PreferenceUpdateReceiver extends BroadcastReceiver {
		@Override 
		public void onReceive(Context context, Intent intent) {
			Log.d(TAG, "Entered onReceive for PreferenceUpdate, Broadcast received");	
			updateActivityPrefs();
		}
		
	}
	
	/*
	 * Task used by the Timer to launch Http requests
	 */
	class TimedHttpTask extends TimerTask {
		@Override
		public void run() {
			Log.d(TAG, "Entered TimedHttpTask, here we need to launch HTTP request");
			// TODO Launch Http task, the function that processes results should launch Json loading task		
		}
	}

}
