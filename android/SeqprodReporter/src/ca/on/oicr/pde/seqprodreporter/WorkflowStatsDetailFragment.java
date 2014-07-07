package ca.on.oicr.pde.seqprodreporter;

import java.text.Format;
import java.text.NumberFormat;
import java.util.LinkedHashMap;

import com.androidplot.ui.YLayoutStyle;
import com.androidplot.ui.YPositionMetric;
import com.androidplot.xy.BarFormatter;
import com.androidplot.xy.BoundaryMode;
import com.androidplot.xy.SimpleXYSeries;
import com.androidplot.xy.ValueMarker.TextOrientation;
import com.androidplot.xy.XValueMarker;
import com.androidplot.xy.XYPlot;
import com.androidplot.xy.YValueMarker;

import android.database.Cursor;
import android.graphics.Paint;
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
	private XYPlot pendingPlot;
	private XYPlot failedPlot;
	String selectedWorkflow;
	
	private static final int BORDER_COLOR = 0xFF000000;
	private static final int HIGHLIGHT_COLOR = 0xFFFFFFFF;
	private static final int FILL_COLOR_COMPLETED = 0XFF000000 ;
	private static final int FILL_COLOR_PENDING = 0xFFCC00FF;
	private static final int FILL_COLOR_FAILED = 0xFFFF0000;
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
			selectedWorkflow = getArguments().getString(ARG_ITEM_ID);
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
		
		completedPlot = (XYPlot) rootView.findViewById(R.id.completePlot);
		completedPlot.setTitle("Completed Workflows");
		
		pendingPlot = (XYPlot) rootView.findViewById(R.id.pendingPlot);
		pendingPlot.setTitle("Pending Workflows");
		
		failedPlot = (XYPlot) rootView.findViewById(R.id.failedPlot);
		failedPlot.setTitle("Failed Workflows");
	/*TODO make function that does it for all 3 XY Plots
		
		for (int workflowType = 0;workflowType<ReporterActivity.types.length;++i){
			setUpXYPlot()
		}
		*/
		int index = 1;
		for (String workflowName : workflowStatsHash.keySet()){
			SimpleXYSeries completedSeries = new SimpleXYSeries(workflowName);
			completedSeries.addFirst(index, workflowStatsHash.get(workflowName)[0]);
			if (workflowName.equals(selectedWorkflow)){
				BarFormatter completedFormat = new BarFormatter(
						HIGHLIGHT_COLOR,BORDER_COLOR);
				completedPlot.addSeries(completedSeries, completedFormat);
			}
			else {
				BarFormatter completedFormat = new BarFormatter(
						FILL_COLOR_COMPLETED,BORDER_COLOR);
				completedPlot.addSeries(completedSeries, completedFormat);
			}
			//TODO: setup marker to indicate either x or Y value
			//XValueMarker valueMarker = new XValueMarker(index,workflowName);
			//valueMarker.setLinePaint(completedPlot.getLayoutManager().getMarginPaint());
			//valueMarker.setTextOrientation(TextOrientation.VERTICAL);
			//valueMarker.setTextPosition(new YPositionMetric(0,YLayoutStyle.ABSOLUTE_FROM_BOTTOM));
			//completedPlot.addMarker(valueMarker);
			++index;
			Log.d(ReporterActivity.TAG, workflowName +": " + workflowStatsHash.get(workflowName)[0]);
		}
		completedPlot.setDomainBoundaries(0, workflowStatsHash.size()+1, BoundaryMode.FIXED);
		completedPlot.setDomainStepValue(1);
		completedPlot.getLegendWidget().setVisible(false);;
		completedPlot.setRangeValueFormat(NumberFormat.getIntegerInstance());
		return rootView;
	}
	
	private void setUpXYPlot(XYPlot plot, int workflowType){
		
	}
}
