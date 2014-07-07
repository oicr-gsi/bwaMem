package ca.on.oicr.pde.seqprodreporter;

import java.util.LinkedHashMap;

import com.androidplot.xy.XYPlot;
import com.example.testingplot.R;

import android.database.Cursor;
import android.os.Bundle;
import android.support.v4.app.Fragment;
import android.util.Log;
import android.view.LayoutInflater;
import android.view.View;
import android.view.ViewGroup;
import android.widget.TextView;
import ca.on.oicr.pde.seqprodprovider.DataContract;
import ca.on.oicr.pde.seqprodreporter.dummy.DummyContent;

/**
 * A fragment representing a single WorkflowStats detail screen. This fragment
 * is either contained in a {@link WorkflowStatsListActivity} in two-pane mode
 * (on tablets) or a {@link WorkflowStatsDetailActivity} on handsets.
 */
public class WorkflowStatsDetailFragment extends Fragment {
	LinkedHashMap<String, Number[]> workflowStatsHash = new LinkedHashMap<String, Number[]>();
	String [] workflowList;
	private XYPlot completedPlot;
	/**
	 * The fragment argument representing the item ID that this fragment
	 * represents.
	 */
	public static final String ARG_ITEM_ID = "item_id";

	/**
	 * The dummy content this fragment is presenting.
	 */

	/**
	 * Mandatory empty constructor for the fragment manager to instantiate the
	 * fragment (e.g. upon screen orientation changes).
	 */
	public WorkflowStatsDetailFragment() {
	}

	@Override
	public void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		if (getArguments().containsKey(ARG_ITEM_ID) && getArguments().containsKey("WorkflowList")) {
			// Load the dummy content specified by the fragment
			// arguments. In a real-world scenario, use a Loader
			// to load content from a content provider.
			workflowList = getArguments().getStringArray("WorkflowList");
			for (int i = 0;i<workflowList.length;++i){
				getWorkflowStats(workflowList[i]);
			}
		}
	}

	private void getWorkflowStats(String workflowName){
		Number[] tmp = new Number[ReporterActivity.types.length];
		for (int i =0;i<ReporterActivity.types.length;++i){
		Cursor c = getActivity().getContentResolver()
			.query(DataContract.CONTENT_URI, new String[]{DataContract.WR_TYPE}, DataContract.WR_TYPE + "=? AND " + DataContract.WORKFLOW + "=? ",new String[]{ReporterActivity.types[i],workflowName} ,null);
		tmp[i] = c.getCount();
		}
		workflowStatsHash.put(workflowName, tmp);
	}
	
	@Override
	public View onCreateView(LayoutInflater inflater, ViewGroup container,
			Bundle savedInstanceState) {
		View rootView = inflater.inflate(
				R.layout.fragment_workflowstats_detail, container, false);
		
		completedPlot = (XYPlot) findViewById(R.id.completePlot);
		completedPlot.setTitle("Completed Workflows");

		return rootView;
	}
}
