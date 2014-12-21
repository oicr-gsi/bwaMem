package ca.on.oicr.pde.seqprodreporter;

import android.app.Fragment;
import android.os.Bundle;
import android.view.LayoutInflater;
import android.view.View;
import android.view.ViewGroup;

public class WorkflowChartFragment extends Fragment {
	private int [] workflowTotals;
	protected int[] getWorkflowTotals() {
		return workflowTotals;
	}

	protected void setWorkflowTotals(int[] workflowTotals) {
		this.workflowTotals = workflowTotals;
	}

	private ChartWidget workflowPieChart; 
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
	public WorkflowChartFragment() {
	}
	
	public WorkflowChartFragment InstanceOf(int[] totals) {
		WorkflowChartFragment fragment = new WorkflowChartFragment();
		fragment.setWorkflowTotals(totals);
		return fragment;
	}

	@Override
	public void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		Bundle arguments = getArguments();
		if (arguments.containsKey(ARG_ITEM_ID)) {
			selectedWorkflow = arguments.getString(ARG_ITEM_ID);
 
			getActivity().setTitle(selectedWorkflow);
			
			if (arguments.containsKey(WorkflowStatsListActivity.WORKFLOW_PIE_CHART_VALUES)){
				workflowTotals = arguments
						.getIntArray(WorkflowStatsListActivity.WORKFLOW_PIE_CHART_VALUES);
			}
			//else if (!selectedWorkflow.equals(WorkflowStatsListFragment.ALL_WORKFLOWS)){
			//	workflowTotals = getPieChartValues();
			//}
		}
	}
	
	
	
	@Override
	public View onCreateView(LayoutInflater inflater, ViewGroup container,
			Bundle savedInstanceState) {
		View rootView = inflater.inflate(
				R.layout.fragment_workflowstats_piechart, container, false);
		
		workflowPieChart = (ChartWidget) rootView.findViewById(R.id.pieChart);
		workflowPieChart.updateTypeData(this.getWorkflowTotals(), null);
		//setUpPieChart();
		return rootView;
	}
	
	public int[] getSelectedWorkflowValues(){
		return workflowTotals;
	}
}
