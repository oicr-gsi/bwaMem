package ca.on.oicr.pde.seqprodreporter;

import java.util.List;

import android.os.Bundle;
import android.support.v4.app.Fragment;
import android.view.LayoutInflater;
import android.view.View;
import android.view.ViewGroup;
import android.widget.ListView;

public class RunListFragment extends Fragment {

	/**
	 * The fragment argument representing the section number for this
	 * fragment.
	 */
	private static final String ARG_SECTION_NUMBER = "section_number";
	private ReportAdapter mAdapter;

	/**
	 * Returns a new instance of this fragment for the given section number.
	 */
	public static RunListFragment newInstance(int sectionNumber) {
		RunListFragment fragment = new RunListFragment();
		Bundle args = new Bundle();
		args.putInt(ARG_SECTION_NUMBER, sectionNumber);
		fragment.setArguments(args);
		return fragment;
	}

	@Override
	public View onCreateView(LayoutInflater inflater, ViewGroup container,
			Bundle savedInstanceState) {
		View rootView = inflater.inflate(R.layout.fragment_reporter,
				container, false);
		ListView listView = (ListView) rootView
				.findViewById(R.id.section_list);
		this.mAdapter = new ReportAdapter(container.getContext());
		int index = this.getSectionNumber() - 1;
		//TODO Need to either load internal file with data or send Http request to get the new json
		new JsonLoaderTask(this, ReporterActivity.types[index]).execute();
		// need to initialize and set the adapter here
		listView.setAdapter(mAdapter);
		return rootView;
	}
	
	public void addLocalReports(List<Report> newReports) {
		this.mAdapter.removeAllViews();
		for (Report r : newReports) {
			this.mAdapter.add(r);
		}
		
		mAdapter.notifyDataSetChanged();
	}
	
	//TODO may need another addRemoteReports function to load data from http
	
	public int getSectionNumber() {
        return getArguments().getInt(ARG_SECTION_NUMBER, 0);
    }

}
