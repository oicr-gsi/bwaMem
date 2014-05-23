package ca.on.oicr.pde.seqprodreporter;

import java.util.Comparator;
import java.util.List;
import java.util.Locale;

import android.os.Bundle;
import android.support.v4.app.Fragment;
import android.text.format.Time;
import android.util.Log;
import android.view.LayoutInflater;
import android.view.View;
import android.view.ViewGroup;
import android.widget.ListView;

public class ReportListFragment extends Fragment {

	/**
	 * The fragment argument representing the section number for this fragment.
	 */
	// private static final String ARG_SECTION_NUMBER = "section_number";
	private ReportAdapter mAdapter;
	private Time lastUpdateTime;
	private int sectionNumber;

	public final Comparator<Report> TIMECOMPARATOR = new ReportTimeComparator();
	public final Comparator<Report> SAMPLECOMPARATOR = new ReportNameComparator();
	public final Comparator<Report> NAMECOMPARATOR = new ReportWorkflowComparator();

	/**
	 * Returns a new instance of this fragment for the given section number.
	 */
	public static ReportListFragment newInstance(int sectionNumber) {
		ReportListFragment fragment = new ReportListFragment();
		fragment.setSectionNumber(sectionNumber);
		return fragment;
	}

	@Override
	public View onCreateView(LayoutInflater inflater, ViewGroup container,
			Bundle savedInstanceState) {

		View rootView = inflater.inflate(R.layout.fragment_reporter, container,
				false);
		ListView listView = (ListView) rootView.findViewById(R.id.section_list);

		int index = this.getSectionNumber() - 1;
		Log.d(ReporterActivity.TAG, "onCreateView called for "
				+ ReporterActivity.types[index]);
		this.mAdapter = new ReportAdapter(container.getContext(), R.layout.fragment_reporter);
		this.mAdapter.setNotifyOnChange(false);
		new JsonLoaderTask(this, ReporterActivity.types[index], this.lastUpdateTime).execute();

		listView.setAdapter(mAdapter);
		return rootView;
	}

	public void addLocalReports(List<Report> newReports) {
		this.mAdapter.removeAllViews();
		//TODO PDE-604 If newReports is empty, add one Report with a TextView text set to "No [type] workflow runs available at this time"
		for (Report r : newReports) {
			this.mAdapter.add(r);
		}
		mAdapter.sortList(TIMECOMPARATOR);
		mAdapter.notifyDataSetChanged();
	}

	// TODO PDE-588 need to switch between these comparators

	public int getSectionNumber() {
		return this.sectionNumber;
	}

	public void setLastUpdateTime(Time t) {
		this.lastUpdateTime = t;
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
			// ascending order
			return Time.compare(repTime1, repTime2);

		}

	};

}
