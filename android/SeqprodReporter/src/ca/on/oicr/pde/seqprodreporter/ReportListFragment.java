package ca.on.oicr.pde.seqprodreporter;

import java.util.Comparator;
import java.util.List;
import java.util.Locale;

import android.app.Activity;
import android.content.SharedPreferences;
import android.os.Bundle;
import android.support.v4.app.Fragment;
import android.text.format.Time;
import android.util.Log;
import android.util.TimeFormatException;
import android.view.LayoutInflater;
import android.view.View;
import android.view.ViewGroup;
import android.widget.ListView;
import android.widget.ProgressBar;

/**
 * A Fragment that contains and shows a list of Reports.
 *
 */
public class ReportListFragment extends Fragment {

	/**
	 * The fragment argument representing the section number for this fragment.
	 */
	private ReportAdapter mAdapter;
	private Time lastUpdateTime;
	private Time firstUpdateTime;
	private int sectionNumber;
	private boolean emptyList;
	private int sortingType;
	private String searchFilter;

	static final int SORT_BY_MODTIME  = 0;
	static final int SORT_BY_WORKFLOW = 1;
	static final int SORT_BY_SAMPLE   = 2;
	// May want to move these comparator's code into a separate class for easy
	// maintenance
	public final Comparator<Report> TIMECOMPARATOR = new ReportTimeComparator();
	public final Comparator<Report> SAMPLECOMPARATOR = new ReportNameComparator();
	public final Comparator<Report> NAMECOMPARATOR = new ReportWorkflowComparator();
    onTimeUpdateListener mCallback;
	/**
	 * Returns a new instance of this fragment for the given section number.
	 * 
	 * @param sectionNumber  is the section number of this fragment
	 * @param searchFilter   search filter String determines which records are getting loaded
	 */
	public static ReportListFragment newInstance(int sectionNumber,
			String searchFilter) {
		ReportListFragment fragment = new ReportListFragment();
		fragment.setSectionNumber(sectionNumber);
		fragment.setSearchFilter(searchFilter);
		return fragment;
	}

	public String getSearchFilter() {
		return searchFilter;
	}
	
	public interface onTimeUpdateListener {
		public void onUpdate(Time updatedTime); 
	}

	/**
	 * This function updates the value of searchFilter variable and updates
	 * 
	 * @param searchFilter (a search query String)
	 */
	public void setSearchFilter(String searchFilter) {
		// THIS CODE RUNS ONLY WHEN ACTIVITY IS VISIBLE
		if (null == this.searchFilter || (!this.searchFilter.equals(searchFilter))) {
			this.searchFilter = searchFilter;
		}
	}

	@Override
	public View onCreateView(LayoutInflater inflater, ViewGroup container,
			Bundle savedInstanceState) {

		View rootView = inflater.inflate(R.layout.fragment_reporter, container,
				false);
		ListView listView = (ListView) rootView.findViewById(R.id.section_list);

		int index = this.getSectionNumber() - 1;
		this.mAdapter = new ReportAdapter(container.getContext(),
				                          R.layout.fragment_reporter);
		this.mAdapter.setNotifyOnChange(false);
		// Restoring last update time, if available:
		this.lastUpdateTime  = new Time();
		this.firstUpdateTime = new Time();
		this.lastUpdateTime.setToNow();
		SharedPreferences sp = getActivity().getSharedPreferences(ReporterActivity.PREFERENCE_FILE, 
				                                                  ReporterActivity.MODE_PRIVATE);
		// need to have just one UpdateTime
		long upTime = sp.getLong("updateLastTime", 0L); // + ReporterActivity.getType(index)
		long crTime = sp.getLong("updateFirstTime" + ReporterActivity.getType(index), 0L);
		String timeRange = sp.getString("pref_summaryScope", "week");
		
		if (upTime > 0) {
			try {
				this.lastUpdateTime.set(upTime);
				if (Time.isEpoch(this.lastUpdateTime)) {
					this.lastUpdateTime.setToNow(); // Prevent passing value created with Time() constructor
				}
			} catch (TimeFormatException tfe) {
				// In case we have a corrupted value, it will be reset in shared preferences
				Log.e(ReporterActivity.TAG, "Time format exception, setting LMT to now...");
				this.lastUpdateTime.setToNow();
			}
		} else {
			this.lastUpdateTime.setToNow();
		}
		
		if (crTime > 0) {
			try {
				this.firstUpdateTime.set(crTime);
				if (Time.isEpoch(this.firstUpdateTime)) {
					this.firstUpdateTime.setToNow(); // Prevent passing value created with Time() constructor
				}
			} catch (TimeFormatException tfe) {
				// In case we have a corrupted value, it will be reset in shared preferences				
				getActivity().getSharedPreferences(ReporterActivity.PREFERENCE_FILE,ReporterActivity.MODE_PRIVATE)
						.edit().putLong("updateFirstTime" + ReporterActivity.getType(this.getSectionNumber() - 1), 0).apply();
			}
		} else {
			this.firstUpdateTime.setToNow();
		}

		new JsonLoaderTask(this, ReporterActivity.getType(index),this.lastUpdateTime).execute(new String[] {getSearchFilter(),timeRange});

		listView.setAdapter(mAdapter);
		return rootView;
	}

	@Override
	public void onAttach(Activity activity) {
		super.onAttach(activity);
		// This makes sure that the container activity has implemented
		// the callback interface. If not, it throws an exception
		try {
			mCallback = (onTimeUpdateListener) activity;
		} catch (ClassCastException e) {
			throw new ClassCastException(activity.toString()
					+ " must implement OnTimeUpdateListener");
		}
	}
	
	@Override
	public void onDestroyView() {
		if (this.lastUpdateTime != null) {
			SharedPreferences sp = getActivity().getSharedPreferences(
					ReporterActivity.PREFERENCE_FILE,
					ReporterActivity.MODE_PRIVATE);

			sp.edit().putLong("updateFirstTime" + ReporterActivity.getType(this.getSectionNumber() - 1),
		            this.firstUpdateTime.toMillis(false)).apply();

		}

		super.onDestroyView();
	}

	/**
	 * This is the function that gets called after JsonLoaderTask finishes
	 * loading data if the received list of reports is empty, a special message
	 * is shown
	 * 
	 * @param newReports
	 */
	public void addLocalReports(List<Report> newReports) {
		ProgressBar a_wheel = (ProgressBar) getView().findViewById(
				R.id.section_progress);
		a_wheel.setVisibility(View.VISIBLE);

		this.mAdapter.removeAllViews();
		if (newReports.size() != 0) {
			for (Report r : newReports) {
				this.mAdapter.add(r);
			}
			this.emptyList = false;
			sortFragment();
		}
		else {
			String type = ReporterActivity.getType(this.sectionNumber-1);
			Report emptyReport = new Report(getString(R.string.empty_message) + " " + type, 
					                        Report.EMPTY_REPORT, "", 0, 0, "", "", "", false);
			this.mAdapter.add(emptyReport);
			this.emptyList = true;
		}
		mAdapter.notifyDataSetChanged();
		a_wheel.setVisibility(View.GONE);
	}

	/**
	 * Invokes the correct comparator to sort a fragment's list based on the
	 * 'sortingType' member variable that is set by the user
	 */
	public void sortFragment() {
		if (this.sortingType == SORT_BY_WORKFLOW) {
			mAdapter.sortList(NAMECOMPARATOR);
		} else if (this.sortingType == SORT_BY_SAMPLE) {
			mAdapter.sortList(SAMPLECOMPARATOR);
		} else {
			mAdapter.sortList(TIMECOMPARATOR);
		}
	}

	public boolean isFragmentListEmpty() {
		return this.emptyList;
	}

	public int getSectionNumber() {
		return this.sectionNumber;
	}

	public ReportAdapter getAdapter() {
		return this.mAdapter;
	}

	public void setSortIndex(int index) {
		this.sortingType = index;
	}

	public void setLastUpdateTime(Time t) {
		if (this.lastUpdateTime == null || this.lastUpdateTime.before(t)) {
			this.lastUpdateTime = t;
			this.mCallback.onUpdate(t);
		}
	}

	public void setFirstUpdateTime(Time t) {
		if (this.firstUpdateTime == null || this.firstUpdateTime.after(t)) {
			this.firstUpdateTime = t;
		}
	}
	
	public Time getFirstUpdateTime() {
		return this.firstUpdateTime;
	}
	
	private void setSectionNumber(int sectionNumber) {
		this.sectionNumber = sectionNumber;
	}

	// Comparators
	private class ReportNameComparator implements Comparator<Report> {
		public int compare(Report report1, Report report2) {

			String repSName1 = report1.getrSampleName().toUpperCase(
					Locale.getDefault());
			String repSName2 = report2.getrSampleName().toUpperCase(
					Locale.getDefault());

			// ascending order
			return repSName1.compareTo(repSName2);

		}
	};

	private class ReportWorkflowComparator implements Comparator<Report> {
		public int compare(Report report1, Report report2) {

			String repWName1 = report1.getrWorkflowName().toUpperCase(
					Locale.getDefault());
			String repWName2 = report2.getrWorkflowName().toUpperCase(
					Locale.getDefault());

			// ascending order
			return repWName1.compareTo(repWName2);

		}

	};

	private class ReportTimeComparator implements Comparator<Report> {
		public int compare(Report report1, Report report2) {

			Time repTime1 = report1.getTimeStamp();
			Time repTime2 = report2.getTimeStamp();

			// newer first
			return Time.compare(repTime2, repTime1);

		}

	};

}
