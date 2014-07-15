package ca.on.oicr.pde.seqprodreporter;

import com.androidplot.pie.PieChart;
import com.androidplot.pie.Segment;
import com.androidplot.pie.SegmentFormatter;
import com.androidplot.ui.SizeLayoutType;
import com.androidplot.ui.SizeMetrics;
import android.database.Cursor;
import android.os.Bundle;
import android.support.v4.app.Fragment;
import android.view.LayoutInflater;
import android.view.View;
import android.view.ViewGroup;
import ca.on.oicr.pde.seqprodprovider.DataContract;
/**
 * A fragment representing a single WorkflowStats detail screen. This fragment
 * is either contained in a {@link WorkflowStatsListActivity} in two-pane mode
 * (on tablets) or a {@link WorkflowStatsDetailActivity} on handsets.
 */
public class WorkflowStatsDetailFragment extends Fragment  {
	
	private int [] workflowTotals;
	private PieChart workflowPieChart; 
	private String selectedWorkflow;	
	/**
	 * The fragment argument representing the item ID that this fragment
	 * represents.
	 */
	public static final String ARG_ITEM_ID = "item_id";
	/**
	 * Mandatory empty constructor for the fragment manager to instantiate the
	 * fragment (e.g. upon screen orientation changes).
	 */
	public WorkflowStatsDetailFragment() {
	}

	@Override
	public void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		Bundle arguments = getArguments();
		if (arguments.containsKey(ARG_ITEM_ID)) {
			selectedWorkflow = arguments.getString(ARG_ITEM_ID);
 
			setActionBarTitle(selectedWorkflow);
			
			if (arguments.containsKey(WorkflowStatsListActivity.WORKFLOW_PIE_CHART_VALUES)){
				workflowTotals = arguments
						.getIntArray(WorkflowStatsListActivity.WORKFLOW_PIE_CHART_VALUES);
			}
			else if (!selectedWorkflow
					.equals(WorkflowStatsListActivity.NO_WORKFLOW_SELECTED)) {
				workflowTotals = getPieChartValues();
			}
		}
	}
	
	@Override
	public View onCreateView(LayoutInflater inflater, ViewGroup container,
			Bundle savedInstanceState) {
		View rootView = inflater.inflate(
				R.layout.fragment_workflowstats_detail, container, false);
		
		workflowPieChart = (PieChart) rootView.findViewById(R.id.pieChart);

		setUpPieChart();

		return rootView;
	}

	private void setActionBarTitle(String title){
		String activityTitle = title
				.equals(WorkflowStatsListActivity.NO_WORKFLOW_SELECTED) 
				? "Workflow Statistics" : title; 
		getActivity().setTitle(activityTitle);
	}
	
	private int[] getPieChartValues(){
		int [] selectedWorkflowNumbers = new int[ReporterActivity.types.length];
			for (int i = 0;i<ReporterActivity.types.length; ++i){
				Cursor c = getActivity().getContentResolver().query(
						DataContract.CONTENT_URI,
						new String[]{DataContract.WORKFLOW},
						DataContract.WORKFLOW + "=? AND " + DataContract.WR_TYPE + "=?",
						new String[]{selectedWorkflow,ReporterActivity.types[i]},
						null);
				
				selectedWorkflowNumbers[i] = c.getCount();
			}
		
		return selectedWorkflowNumbers;
	}
	
	private void setUpPieChart(){
		workflowPieChart.getTitleWidget().getLabelPaint().setTextSize(40f);
		workflowPieChart.getTitleWidget().setSize(
				new SizeMetrics(100f, SizeLayoutType.ABSOLUTE,
						150f, SizeLayoutType.FILL));
	
		for (int i = 0;i<ReporterActivity.types.length;++i){
			if (0 != workflowTotals[i]){
				String segmentLabel = ReporterActivity.types[i].substring(0, 1).toUpperCase() 
						+ ReporterActivity.types[i].substring(1) + ": " + workflowTotals[i];
				
				Segment pieSegment = new Segment(segmentLabel,(Number) workflowTotals[i]);
				SegmentFormatter segmentFormatter = new SegmentFormatter(getTypeFillColor(i), 
						getResources().getColor(R.color.black));
				
				segmentFormatter.getLabelPaint().setColor(getResources().getColor(R.color.black));
				segmentFormatter.getLabelPaint().setTextSize(20f);
				workflowPieChart.addSegment(pieSegment, segmentFormatter);
			}
		}
	}
	
	public int getTypeFillColor(int workflowType){
		int fillColor;
		switch(workflowType){
		case 0:
			fillColor = getResources().getColor(R.color.green);
			break;
		case 1:
			fillColor = getResources().getColor(R.color.red);
			break;
		case 2: 
			fillColor = getResources().getColor(R.color.yellow);
			break;
		default:
			fillColor = getResources().getColor(R.color.black);
			break;
	}
		return fillColor;
	}

	
}
