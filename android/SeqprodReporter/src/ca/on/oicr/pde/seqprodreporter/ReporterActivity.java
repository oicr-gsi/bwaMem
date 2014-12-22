package ca.on.oicr.pde.seqprodreporter;

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Timer;
import java.util.TimerTask;

import android.annotation.SuppressLint;
import android.app.AlertDialog;
import android.app.Notification;
import android.app.NotificationManager;
import android.app.PendingIntent;
import android.app.SearchManager;
import android.content.BroadcastReceiver;
import android.content.Context;
import android.content.DialogInterface;
import android.content.Intent;
import android.content.IntentFilter;
import android.content.SharedPreferences;
import android.database.Cursor;
import android.graphics.Color;
import android.media.RingtoneManager;
import android.os.Bundle;
import android.support.v4.app.Fragment;
import android.support.v4.app.FragmentManager;
import android.support.v4.app.FragmentPagerAdapter;
import android.support.v4.app.FragmentTransaction;
import android.support.v4.content.LocalBroadcastManager;
import android.support.v4.view.MenuItemCompat;
import android.support.v4.view.ViewPager;
import android.support.v7.app.ActionBar;
import android.support.v7.app.ActionBarActivity;
import android.support.v7.widget.SearchView.OnCloseListener;
import android.text.format.Time;
import android.util.Log;
import android.util.TimeFormatException;
import android.view.Menu;
import android.view.MenuItem;
import android.view.View;
import android.view.ViewGroup.OnHierarchyChangeListener;
import android.widget.TextView;
import android.widget.Toast;
import ca.on.oicr.pde.seqprodprovider.DataContract;

/**
 * The launcher Activity for the App. This Activity contains three
 * ReportListFragments as well as a ViewPager that allows the user to navigate
 * between the three ReportListFragments. This extended ActionBarActivity class
 * also implements a ActionBar.TabListener, which allows the user to see the
 * corresponding ReportListFragment on select of a tab. Three tabs exist at the
 * moment ("completed", "failed", "pending"), the TabListener is implemented in
 * such a way so that on select of a tab the user will be navigated to the
 * ReportListFragment that contains the list of workflows that are of the type
 * represented by the selected tab.
 *
 * @see ActionBarActivity
 * @see ViewPager
 */
public class ReporterActivity extends ActionBarActivity implements
		ActionBar.TabListener, ReportListFragment.onTimeUpdateListener {

	/**
	 * The {@link android.support.v4.view.PagerAdapter} that will provide
	 * fragments for each of the sections. We use a {@link FragmentPagerAdapter}
	 * derivative, which will keep every loaded fragment in memory. If this
	 * becomes too memory intensive, it may be best to switch to a
	 * {@link android.support.v4.app.FragmentStatePagerAdapter}.
	 */
	SectionsPagerAdapter mSectionsPagerAdapter;
	private final PreferenceUpdateReceiver prefUpdateReceiver = new PreferenceUpdateReceiver();
	private final DataUpdateReceiver dataUpdateReceiver = new DataUpdateReceiver();
	/**
	 * The {@link ViewPager} that will host the section contents.
	 */
	ViewPager mViewPager;

	public static final String TAG = "Seqprodbio Reporter";

	protected static String[] types = { "completed", "failed", "pending" };
	protected static String SYNC_OFF;
	protected static final String TIMER_NAME = "Seqprodbio Timer";
	protected static final String PREFERENCE_FILE = "seqprod.conf";
	private static final long INITIAL_TIMER_DELAY = 5 * 1000L;

	private static final String NOTIFICATIONS_OFF = "Off";
	private static final String NOTIFICATIONS_WEB_UPDATES = "Web Updates";
	private static final String NOTIFICATIONS_CRITICAL_UPDATES = "Critical Updates";
	private static final String NOTIFICATIONS_CRITICAL_UPDATES_SOUND = "Critical Updates With Sound";
	protected static final int COMPLETED_WORKFLOW_TAB_INDEX = 0;
	protected static final int FAILED_WORKFLOW_TAB_INDEX = 1;
	protected static final int PENDING_WORKFLOW_TAB_INDEX = 2;

	static final String PREFCHANGE_INTENT = "ca.on.oicr.pde.seqprodreporter.prefsChanged";
	static final String DATACHANGE_INTENT = "ca.on.oicr.pde.seqprodreporter.updateLoaded";
	protected static final long [] updateRanges = {7 * 24 * 3600 * 1000L,		//WEEK
		                                           30 * 24 * 3600 * 1000L,		//MONTH
                                                   365 * 24 * 3600 * 1000L,		//YEAR
                                                   10 * 365 * 24 * 3600 * 1000L};//DECADE

	private String updateHost;
	private String updateRange;
	private String notificationSetting;
	private String mSearchQuery;
	private boolean timerScheduled = false;
	private boolean isVisible;
	private boolean[] dataRangeRequested;
	private boolean refreshEnabled = true;
	private Time lastModifiedFailedTime;
	private Time lastModifiedCachedTime;
	private int sortIndex;
	private int mCurrentTabIndex;
	private int updateFrequency; // in minutes

	private SharedPreferences sp;
	private Timer timer;

	@Override
	protected void onCreate(Bundle savedInstanceState) {

		super.onCreate(savedInstanceState);
		SYNC_OFF = getResources().getString(
				R.string.pref_automaticUpdates_default);

		setContentView(R.layout.activity_reporter);

		((MainApplication) getApplication()).setisCurrentActivityVisible(true);

		// Register receivers for preference and data updates
		LocalBroadcastManager lmb = LocalBroadcastManager.getInstance(this);
		IntentFilter prefchangeFilter = new IntentFilter(PREFCHANGE_INTENT);
		lmb.registerReceiver(prefUpdateReceiver, prefchangeFilter);
		IntentFilter datachangeFilter = new IntentFilter(DATACHANGE_INTENT);
		lmb.registerReceiver(dataUpdateReceiver, datachangeFilter);

		this.sp = getSharedPreferences(PREFERENCE_FILE, MODE_PRIVATE);
		// Read preferences
		this.updateActivityPrefs();
		// Set up the action bar.
		final ActionBar actionBar = getSupportActionBar();
		actionBar.setNavigationMode(ActionBar.NAVIGATION_MODE_TABS);
		actionBar.setDisplayShowHomeEnabled(true);
		actionBar.setDisplayShowTitleEnabled(true);

		// Create the adapter that will return a fragment for each of the three
		// primary sections of the activity.
		mSectionsPagerAdapter = new SectionsPagerAdapter(
				getSupportFragmentManager());
		// Set up the ViewPager with the sections adapter.
		mViewPager = (ViewPager) findViewById(R.id.pager);
		mViewPager.setAdapter(mSectionsPagerAdapter);

		// Allows the 2 other tab's fragments that are in idle state to be
		// loaded
		// alongside the current selected tab's fragment
		mViewPager.setOffscreenPageLimit(types.length - 1);

		/*
		 * When swiping between different sections, select the corresponding
		 * tab. We can also use ActionBar.Tab#select() to do this if we have a
		 * reference to the Tab.
		 */
		mViewPager
				.setOnPageChangeListener(new ViewPager.SimpleOnPageChangeListener() {
					@Override
					public void onPageSelected(int position) {
						actionBar.setSelectedNavigationItem(position);
					}
				});

		/*
		 * This is needed for selecting proper tab when restoring activity or
		 * returning from Stats activity - it takes care of selecting tab as it
		 * gets added (May be done differently, perhaps)
		 */
		mViewPager
				.setOnHierarchyChangeListener(new OnHierarchyChangeListener() {

					@Override
					public void onChildViewAdded(View parent, View child) {
						if (mCurrentTabIndex != 0
								&& mViewPager.getChildCount() >= mCurrentTabIndex
								&& mCurrentTabIndex != mViewPager
										.getCurrentItem())
							mViewPager.setCurrentItem(mCurrentTabIndex);
					}

					@Override
					public void onChildViewRemoved(View parent, View child) {
						// Do nothing, this is not supposed to happen
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

		if (null == lastModifiedFailedTime)
			lastModifiedFailedTime = new Time();
		if (null == lastModifiedCachedTime)
			lastModifiedCachedTime = new Time();

	}

	@Override
	protected void onPostCreate(Bundle savedInstanceState) {
		super.onPostCreate(savedInstanceState);
		// Handle intent that comes from Statistics Activity
		Intent startIntent = this.getIntent();
		if (null != startIntent)
			handleIntent(startIntent);
	}

	@Override
	protected void onPause() {
		((MainApplication) getApplication()).setisCurrentActivityVisible(false);
		this.isVisible = false;
		// store cachedLastModTime in sp
		Log.d(TAG,"Setting on pause, overwriting LMT in SharedPreferences with cached value...");
		sp.edit().putLong("updateLastTime", 
				          this.lastModifiedCachedTime.toMillis(false)).apply();
		super.onPause();
	}

	@Override
	protected void onResume() {
		((MainApplication) getApplication()).setisCurrentActivityVisible(true);
		// restore time stored in 
		updateUiLMT(sp.getString("updateTime", ""));
		this.isVisible = true;
		this.lastModifiedCachedTime.set(sp.getLong("updateLastTime", 0L));
		// update fragments when searching for a query (Device orientation
		// change handled elsewhere)
		if (!mSectionsPagerAdapter.fragments.isEmpty()) {
			List<ReportListFragment> fragments = mSectionsPagerAdapter.fragments;
			for (int i = 0; i < fragments.size(); i++) {
				ReportListFragment tmp = fragments.get(i);
				if (null != tmp)
					tmp.setSearchFilter(this.mSearchQuery);
			}

			mSectionsPagerAdapter.notifyDataSetChanged();
		}
		super.onResume();
	}

	@Override
	protected void onSaveInstanceState(Bundle outState) {
		Log.d(TAG, "Saving Instance...");
		if (null != this.lastModifiedFailedTime)
			outState.putString("lastModifiedFailedTime",
					lastModifiedFailedTime.format2445());
		if (null != this.dataRangeRequested)
			outState.putBooleanArray("dataRangeRequested",
					this.dataRangeRequested);
		outState.putInt("currentlySelectedTab", this.mCurrentTabIndex);
		if (null != this.mSearchQuery && !this.mSearchQuery.isEmpty())
			outState.putString("currentSearchQuery", this.mSearchQuery);
	}

	@Override
	protected void onRestoreInstanceState(Bundle savedInstanceState) {
		Log.d(TAG, "Restoring Instance...");
		super.onRestoreInstanceState(savedInstanceState);
		if (null != savedInstanceState) {
			try {
				this.lastModifiedFailedTime = new Time();
				this.lastModifiedFailedTime.parse(savedInstanceState
						.getString("lastModifiedFailedTime"));
				this.mCurrentTabIndex = savedInstanceState
						.getInt("currentlySelectedTab");
				this.mSearchQuery = savedInstanceState
						.getString("currentSearchQuery");
				this.dataRangeRequested = savedInstanceState
						.getBooleanArray("dataRangeRequested");
			} catch (Exception e) {
				Log.d(TAG, "Last Failed Time could not be retrieved");
			}
		} else {
			this.mCurrentTabIndex = 0;
		}
	}

	@Override
	public boolean onCreateOptionsMenu(Menu menu) {
		// Inflate the menu; this adds items to the action bar if it is present.
		getMenuInflater().inflate(R.menu.reporter, menu);
		MenuItem searchItem = menu.findItem(R.id.action_search);
		SeqprodSearchView searchView = (SeqprodSearchView) MenuItemCompat
				.getActionView(searchItem);
		if (null != searchView) {
			SearchManager searchManager = (SearchManager) getSystemService(Context.SEARCH_SERVICE);
			// Assumes current activity is the search-able activity
			searchView.setSearchableInfo(searchManager
					.getSearchableInfo(getComponentName()));
			searchView.setIconifiedByDefault(true);
			searchView.setSearchString(this.mSearchQuery);
			if (this.mSearchQuery != null && !this.mSearchQuery.isEmpty()) {
				searchView.setQuery(this.mSearchQuery, false);
			}

			searchView.setOnCloseListener(new OnCloseListener() {
				@Override
				public boolean onClose() {
					// Indicate that we handled this request and don't want to
					// reset - helpful when we want to get to the tabs
					return mSearchQuery.isEmpty() ? false : true;
				}
			});

		}
		return super.onCreateOptionsMenu(menu);
	}

	private boolean isUpdateHostSet() {
		return !(null == this.updateHost || this.updateHost.isEmpty());
	}

	private boolean isUpdateRangeSet() {
		return !(null == this.updateRange || this.updateRange.isEmpty());
	}

	@Override
	public boolean onOptionsItemSelected(MenuItem item) {
		int id = item.getItemId();
		if (id == R.id.action_settings) {
			Intent setPrefs = new Intent(this, SeqprodPreferencesActivity.class);
			startActivity(setPrefs);
			return true;
		}
		// An alert dialog is invoked for the user to select the sorting method
		// to apply on the list of each fragment
		else if (id == R.id.sort_by) {
			AlertDialog.Builder builder = new AlertDialog.Builder(this);
			builder.setTitle(R.string.sort_dialog).setSingleChoiceItems(
					R.array.sorting_method_types, sortIndex,
					new DialogInterface.OnClickListener() {
						public void onClick(DialogInterface dialog, int selected) {
							List<ReportListFragment> fragments = mSectionsPagerAdapter.fragments;
							sortIndex = selected;
							for (int i = 0; i < fragments.size(); ++i) {
								ReportListFragment tmp = fragments.get(i);
								tmp.setSortIndex(selected);
								tmp.sortFragment();
								tmp.getAdapter().notifyDataSetChanged();
							}
							dialog.dismiss();
						}
					});
			builder.show();
		} else if (id == R.id.action_refresh) {
			if (!this.refreshEnabled) {
				Log.d(TAG, "Refresh disabled, won't send request...");
				return true;
			}
			this.refreshEnabled = false;
			// Instance where a timer is set with automatic updates
			// Timer is cancelled, and re-scheduled with the same interval but
			// starts the instance
			// when the icon is pressed with no delay
			if (null != timer) {
				Toast.makeText(ReporterActivity.this,
						"Lists Are Being Refreshed", Toast.LENGTH_LONG).show();
				long INTERVAL = this.updateFrequency * 60 * 1000L;
				this.timer.cancel();
				this.timer = new Timer();
				// record LUT here (regardless of visibility)
				this.lastModifiedCachedTime.set(sp.getLong("updateLastTime", 0L));
				this.timer.schedule(new TimedHttpTask(), 0, INTERVAL);
			} else {
				// Instance where if either updateHost and updateRange are not
				// set in the preferences
				// an alert dialog will appear telling the user to set the
				// corresponding fields
				if (!isUpdateHostSet() || !isUpdateRangeSet()) {
					AlertDialog.Builder builder = new AlertDialog.Builder(this);
					builder.setTitle(R.string.refresh_error_title).setMessage(
							R.string.refresh_error_message);
					builder.setPositiveButton("Ok",
							new DialogInterface.OnClickListener() {
								public void onClick(DialogInterface dialog,
										int select) {
									dialog.dismiss();
								}
							}).show();
				}
				// Instance where the preferences are filled in, however the
				// automatic updates are turned off
				// Will just schedule a single execution of TimedHttpTask with
				// no change in the preferences
				else {
					Toast.makeText(ReporterActivity.this,
							"Lists Are Being Refreshed", Toast.LENGTH_LONG)
							.show();
					new Timer().schedule(new TimedHttpTask(), 0);
				}
			}
		} else if (id == R.id.stats) {
			if (isDatabaseEmpty()) {
				AlertDialog.Builder builder = new AlertDialog.Builder(this);
				builder.setTitle(R.string.stats_error_title).setMessage(
						R.string.stats_error_message);
				builder.setPositiveButton("Ok",
						new DialogInterface.OnClickListener() {
							public void onClick(DialogInterface dialog,
									int select) {
								dialog.dismiss();
							}
						}).show();
			} else {
				Intent intent = new Intent(this,
						WorkflowStatsActivity.class);
						//WorkflowStatsListActivity.class);
				startActivity(intent);
				return true;
			}
		}
		return super.onOptionsItemSelected(item);
	}

	// Checks if the database is empty by using cursor results, used when we
	// shouldn't be able to go to stats page
	private boolean isDatabaseEmpty() {
		Cursor c = this.getContentResolver().query(DataContract.CONTENT_URI,
				new String[] { DataContract.WORKFLOW }, null, null, null);
		boolean isEmpty = !c.moveToFirst();
		c.close();
		return isEmpty;
	}

	@Override
	public void onTabSelected(ActionBar.Tab tab,
			FragmentTransaction fragmentTransaction) {
		// When the given tab is selected, switch to the corresponding page in
		// the ViewPager.
		mViewPager.setCurrentItem(tab.getPosition());
		mCurrentTabIndex = tab.getPosition();
	}

	@Override
	public void onTabUnselected(ActionBar.Tab tab,
			FragmentTransaction fragmentTransaction) {
	}

	@Override
	public void onTabReselected(ActionBar.Tab tab,
			FragmentTransaction fragmentTransaction) {
	}

	@Override
	protected void onNewIntent(Intent intent) {
		setIntent(intent);
		handleIntent(intent);
	}

	/**
	 * handleIntent
	 * <p>
	 * Introduced to handle SearchView widget events and Intent from Stats
	 * Activity
	 * <p>
	 * Gets all fragments and set their search Filter to query, entered in
	 * SearchView
	 * 
	 * @param Intent
	 */

	private void handleIntent(Intent intent) {
		// Search widget intent
		if (Intent.ACTION_SEARCH.equals(intent.getAction())) {
			this.isVisible = true;
			Log.d(TAG, "Calling Search/Filter code...");
			this.mSearchQuery = intent.getStringExtra(SearchManager.QUERY);
			// Statistics Activity event
		} else if (Intent.ACTION_RUN.equals(intent.getAction())) {
			Object starter = intent.getExtras().get("selectedTab");
			if (null != starter) {
				try {
					int selectedTab = Integer.parseInt(starter.toString());
					if (selectedTab >= 0 && selectedTab < types.length) {
						this.mCurrentTabIndex = selectedTab;
					}
				} catch (NumberFormatException nfe) {
					Log.e(TAG, "Received invalid tab index from Intent");
				}
			}
		} else { // Debugging only
			String act = intent.getAction();
			Log.d(TAG, "Action passed: " + act);
		}
	}

	public static String getType(int index) {
		return types[index];
	}

	/**
	 * A function that returns a formatted user-friendly Time representation.
	 * 
	 * @param Time
	 *            time
	 */
	public static String timeToStringConverter(Time time) {
		return time.format("%Y-%m-%d %H:%M:%S");
	}
	
	
	/**
	 *  A function for getting earliest record's time in milliseconds
	 *  for a specified time range (week, month etc)
	 */
	protected static long getEarliestMillis(String timeRange) {

		long earliest =  Long.valueOf(System.currentTimeMillis());	
		if (null != timeRange && !timeRange.isEmpty()) {	    
		    // Not clear at this point what is the best way to handle this other than
		    // by using long if-else block
		    if (timeRange.equals("decade")) {
		      earliest -= updateRanges[3];		    	
		    } else if (timeRange.equals("year")) {
		      earliest -= updateRanges[2];
		    } else if (timeRange.equals("month")) {
		      earliest -= updateRanges[1];
		    } else { // default to 'week'
		      earliest -= updateRanges[0];
		    }	    
		} else {
			earliest -= updateRanges[0];
		}
		
		return earliest;
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
			if (position >= this.fragments.size() || this.fragments.size() == 0
					|| null == this.fragments.get(position)) {

				fragments.add(position, ReportListFragment.newInstance(
						position + 1, mSearchQuery));
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

		// Forces each Item in the list to be re-created which allows the lists
		// to be updated dynamically

		@Override
		public int getItemPosition(Object object) {
			return POSITION_NONE;
		}

	}

	/*
	 * This function sets the Timer, but DOES NOT launches Http task
	 */
	private void scheduleUpdate() {
		// Update data range here as well
		Log.v(TAG, "Entered sheduleUpdates");
		// Start update ONLY if host URL is valid
		if (!isUpdateHostSet() || !isUpdateRangeSet()
				|| this.updateFrequency == 0) {
			if (null != this.timer) {
				this.timer.cancel();
				this.timer = null;
				Log.d(TAG, "Http Updates Canceled");
			}
			this.timerScheduled = false;
			this.updateFrequency = 0;
			return;
		}
		Log.d(TAG, "Will schedule Timer to trigger updates from "
				+ this.updateHost);
		// Schedule Alert here
		long INTERVAL = this.updateFrequency * 60 * 1000L;
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
		String uf = sp.getString("pref_syncFreq", SYNC_OFF);
		if (!uf.equals(SYNC_OFF)) {
			try {
				this.updateFrequency = Integer.parseInt(uf.substring(0,
						uf.lastIndexOf(" ")));
			} catch (NumberFormatException ne) {
				Log.e(TAG,
						"Malformed value for update frequency, will not update");
				this.updateFrequency = 0;
			}
		} else {
			this.updateFrequency = 0;
		}

		// PDE-650 handle time ranges other than 'week' here
		this.updateRange = sp.getString("pref_summaryScope", null);
		if (null == this.updateRange || null == this.updateHost)
			return;

		if (!this.updateRange.equals(getResources().getStringArray(
				R.array.pref_summaryScope_entries)[0])
				&& this.updateFrequency != 0)
			requestLargeUpdate();

		this.notificationSetting = sp.getString("pref_notificationSettings",
				NOTIFICATIONS_WEB_UPDATES);

		scheduleUpdate();
	}

	/**
	 * Update Last Modified Time in UI
	 * 
	 * @param updateTime
	 */
	private void updateUiLMT(String updateTime) {
		if (!updateTime.equals("")) {
			TextView updateView = (TextView) findViewById(R.id.updateTimeView);
			updateView.setVisibility(View.VISIBLE);
			Log.d(TAG, "Update Time set to " + updateTime);
			String newTitle = "Most Recent Workflow Modification Time: "
					+ updateTime;
			updateView.setText(newTitle);
			updateView.invalidate();
		}
	}

	/**
	 * Function for sending HTTP requests for data range other than 'week'
	 */
	private void requestLargeUpdate() {
		// check dataRangeRequested flag
		String[] ranges = getResources().getStringArray(
				R.array.pref_summaryScope_entries);
		int rangeIndex = 0;
		for (int i = 1; i < ranges.length; i++) { // we don't check 'week' which
													// is the element with index
													// 0
			if (this.updateRange.equals(ranges[i]))
				rangeIndex = i;
		}
		if (rangeIndex == 0)
			return;
		
		Log.d(TAG, "Update range longer than a week chosen");
		
		if (null == this.dataRangeRequested)
			this.dataRangeRequested = new boolean[ranges.length];
		else if (this.dataRangeRequested[rangeIndex])
			return;

		// check the oldest (completed) fragments' firstUpdated value
		Log.d(TAG, "Checking if we already have data in requested range...");
		long rangeLimit = getEarliestMillis(this.updateRange); 

		try {
			ReportListFragment f = (ReportListFragment) mSectionsPagerAdapter
					.getItem(0);
			if (null == f.getFirstUpdateTime()
					|| f.getFirstUpdateTime().toMillis(false) <= rangeLimit) {
				this.dataRangeRequested[rangeIndex] = true;
				return;
			}
			// if all checks negative, launch http request
			Log.d(TAG, "Requesting a large update...");
			new getreportHTTP(getApplicationContext(), this.updateHost,
					this.updateRange).execute();
			// update dataRangeRequested
			this.dataRangeRequested[rangeIndex] = true;
		} catch (NullPointerException npe) {
			Log.e(TAG, "Error getting first update time, won't request data");
		}
	}

	/*
	 * Broadcast Receiver for Preference Update Broadcast, updates variables
	 * with preference values
	 */
	class PreferenceUpdateReceiver extends BroadcastReceiver {
		@Override
		public void onReceive(Context context, Intent intent) {
			Log.d(TAG, "Entered onReceive for PreferenceUpdateReceiver");
			updateActivityPrefs();
		}

	}

	@SuppressLint("NewApi")
	/*
	 * Broadcast Receiver for Data Update Broadcast
	 */
	class DataUpdateReceiver extends BroadcastReceiver {
		@Override
		public void onReceive(Context context, Intent intent) {
			if (isVisible != ((MainApplication) getApplication())
					.getisCurrentActivityVisible()) {
				LocalBroadcastManager.getInstance(getApplicationContext())
						.unregisterReceiver(this);
			}

			else {
				// The intent will only have the extra "updateTime" when the
				// reports have been successfully loaded from the
				// GetReportHttp async task, therefore this will notify the user
				// that there was an error loading the reports
				// depending on if they are on screen (Toast message) or off
				// screen (Notification)
				if (!intent.hasExtra("updateTime")) {
					if (ReporterActivity.this.isVisible) {
						Toast.makeText(ReporterActivity.this,
								"Error: No New Reports Could be Loaded",
								Toast.LENGTH_LONG).show();
					} else if (!notificationSetting.equals(NOTIFICATIONS_OFF)) {
						Notification.Builder notificationBuilder = setUpNotificationBuilder(context);
						notificationBuilder
								.setTicker("Error: No New Reports Could be Loaded")
								.setContentText("An Error Occured While Loading The Reports");
						passNotificationBuildertoManager(context,
								notificationBuilder);
					}
				}

				else {
					sp.edit()
							.putString("updateTime", intent.getStringExtra("updateTime"))
							.apply();
					Log.d(TAG, "Entered onReceive for DataUpdate, Broadcast received");

					if (ReporterActivity.this.isVisible) {
						Toast.makeText(ReporterActivity.this,
								"Update Received", Toast.LENGTH_SHORT).show();
						updateUiLMT(intent.getStringExtra("updateTime"));
					} else {
						if (!notificationSetting.equals(NOTIFICATIONS_OFF)) {
							Notification.Builder notificationBuilder = setUpNotificationBuilder(context);

							if (notificationSetting
									.equals(NOTIFICATIONS_WEB_UPDATES)) {
								notificationBuilder
										.setTicker("Update Received")
										.setContentText("Update Received");

							} else if (!Time.isEpoch(lastModifiedFailedTime)
									&& isFailedModified(intent)
									&& (notificationSetting
											.equals(NOTIFICATIONS_CRITICAL_UPDATES) || notificationSetting
											.equals(NOTIFICATIONS_CRITICAL_UPDATES_SOUND))) {
								notificationBuilder
										.setTicker("Critical Update Received")
										.setContentText(
												"Critical Update Received")
										.setLights(Color.RED, 500, 1000);
								if (notificationSetting
										.equals(NOTIFICATIONS_CRITICAL_UPDATES_SOUND)) {
									// May need to change the notification sound
									// type
									notificationBuilder
											.setSound(RingtoneManager
													.getDefaultUri(RingtoneManager.TYPE_NOTIFICATION));
								}
							}
							// Pass the Notification to the NotificationManager:
							passNotificationBuildertoManager(context,
									notificationBuilder);
						}
					}

					if (isFailedModified(intent)) {
						lastModifiedFailedTime.parse(intent
								.getStringExtra("modifiedFailedTime"));
					}

					if (ReporterActivity.this.isVisible)
						mSectionsPagerAdapter.notifyDataSetChanged();

				}
			}
			ReporterActivity.this.refreshEnabled = true;
		}

		private Notification.Builder setUpNotificationBuilder(Context context) {
			Intent mNIntent = new Intent(context, ReporterActivity.class);
			PendingIntent mCIntent = PendingIntent.getActivity(context, 0,
					mNIntent, Intent.FLAG_ACTIVITY_NEW_TASK);
			Notification.Builder notificationBuilder = new Notification.Builder(
					context);
			notificationBuilder
					.setSmallIcon(android.R.drawable.stat_sys_warning)
					.setAutoCancel(true).setContentTitle("Seqprod Reporter")
					.setContentIntent(mCIntent);
			return notificationBuilder;
		}

		private void passNotificationBuildertoManager(Context context,
				Notification.Builder notificationBuilder) {
			NotificationManager mNotificationManager = (NotificationManager) context
					.getSystemService(Context.NOTIFICATION_SERVICE);
			mNotificationManager.notify(0, notificationBuilder.build());
		}

		private boolean isFailedModified(Intent intent) {
			return intent.hasExtra("modifiedFailedTime");
		}
	}

	/**
	 * Class used by the Timer to launch Http requests
	 */
	class TimedHttpTask extends TimerTask {
		@Override
		public void run() {
			if (isVisible != ((MainApplication) getApplication())
					.getisCurrentActivityVisible()) {
				timer.cancel();
			} else {
				Log.d(TAG,
						"Entered TimedHttpTask, about to launch HTTP request, caching LMT...");
				//Cache update time before launching HTTP update
				long storedMillis = sp.getLong("updateLastTime", 0L);
				
				if (storedMillis > 0)
					lastModifiedCachedTime.set(sp.getLong("updateLastTime", 0L));
				// We ask for 'week' update here and launch getreportHTTP
				new getreportHTTP(getApplicationContext(), sp.getString(
						"pref_hostName", null), getResources().getStringArray(
						R.array.pref_summaryScope_entries)[0])
						.execute(lastModifiedFailedTime);
			}
		}

	}

	@Override
	public synchronized void onUpdate(Time updatedTime) {

		if (!((MainApplication) getApplication()).getisCurrentActivityVisible())
			return;

		Log.d(TAG, "LMT is getting checked against the stored value...");
		// this function will update sp-stored time if it is after the
		// stored value DO NOT touch cached time here!
		long upTime = sp.getLong("updateLastTime", 0L);
		Time restoredTime = new Time();
		try {
			restoredTime.set(upTime);
			if (updatedTime.after(restoredTime)) {
				Log.d(TAG, "Updating...");
				this.sp.edit()
					   .putLong("updateLastTime", updatedTime.toMillis(false))
					   .apply();
			}
		} catch (TimeFormatException tfe) {
			Log.d(TAG, "LMT wasn't updated b/c the value is wrong...");
		}

	}

}
