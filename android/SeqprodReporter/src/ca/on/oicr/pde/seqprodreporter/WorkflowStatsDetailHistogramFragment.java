package ca.on.oicr.pde.seqprodreporter;

import java.text.NumberFormat;
import java.util.LinkedHashMap;
import java.util.List;

import com.androidplot.ui.YLayoutStyle;
import com.androidplot.ui.YPositionMetric;
import com.androidplot.xy.BarFormatter;
import com.androidplot.xy.BarRenderer;
import com.androidplot.xy.BoundaryMode;
import com.androidplot.xy.SimpleXYSeries;
import com.androidplot.xy.XValueMarker;
import com.androidplot.xy.XYPlot;
import com.androidplot.xy.XYSeriesRenderer;

import android.support.v4.app.Fragment;
import android.os.Bundle;
import android.view.LayoutInflater;
import android.view.View;
import android.view.ViewGroup;

public class WorkflowStatsDetailHistogramFragment extends Fragment {
	private XYPlot completedPlot;
	private XYPlot pendingPlot;
	private XYPlot failedPlot;
	private LinkedHashMap<String, int[]> workflowStatsHash;
	private String selectedWorkflow;
	private WorkflowStatsListActivity activity;
	
	private static final double VALUE_LABEL_OFFSET = 0.1;
	private static final int HIGHLIGHT_COLOR = 0xFFFFFFFF;
	private static final float GRAPH_TOP_MARGIN = 15f;
	private static final float BAR_WIDTH = 40f;
	
	
	public WorkflowStatsDetailHistogramFragment(){
	}
	
	@Override
	public void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		
		activity = (WorkflowStatsListActivity) getActivity();
		workflowStatsHash = activity.getWorkflowStatsHash();
		
		if (getArguments().containsKey(WorkflowStatsDetailFragment.ARG_ITEM_ID)) {
			selectedWorkflow = getArguments().getString(WorkflowStatsDetailFragment.ARG_ITEM_ID);
		}
	}
	
	@Override
	public View onCreateView(LayoutInflater inflater, ViewGroup container,
			Bundle savedInstanceState) {
		View rootView = inflater.inflate(
				R.layout.fragment_workflowstats_histogramdetail, container, false);
		
		completedPlot = (XYPlot) rootView.findViewById(R.id.completePlot);
		completedPlot.setTitle("Completed Workflows");
		
		failedPlot = (XYPlot) rootView.findViewById(R.id.failedPlot);
		failedPlot.setTitle("Failed Workflows");
		
		pendingPlot = (XYPlot) rootView.findViewById(R.id.pendingPlot);
		pendingPlot.setTitle("Pending Workflows");

		setUpXYPlot(completedPlot, 0);
		setUpXYPlot(failedPlot, 1);
		setUpXYPlot(pendingPlot, 2);

		return rootView;
	}
	
	private void setUpXYPlot(XYPlot plot, int workflowType){
		int workflowIndex = 1;
		int fillColor = WorkflowStatsDetailFragment.getTypeFillColor(workflowType);
		int maxRange = 0;
		
		for (String workflowName : workflowStatsHash.keySet()){

			if (1 == workflowIndex){
				maxRange = workflowStatsHash.get(workflowName)[workflowType];
			}
			else if (maxRange < workflowStatsHash.get(workflowName)[workflowType]){
				maxRange = workflowStatsHash.get(workflowName)[workflowType];
			}
			
			BarFormatter barFormatter;
			SimpleXYSeries series = new SimpleXYSeries(workflowName);
			series.addFirst(workflowIndex, 
					workflowStatsHash.get(workflowName)[workflowType]);
			
			if (!selectedWorkflow.equals("NoSelectedWorkflow") 
					&& workflowName.equals(selectedWorkflow)){
				plot.addMarker(new XValueMarker(workflowIndex- VALUE_LABEL_OFFSET,
						((Number)workflowStatsHash.get(workflowName)[workflowType]).toString(),
							new YPositionMetric(0,YLayoutStyle.ABSOLUTE_FROM_BOTTOM),
								0,WorkflowStatsDetailFragment.BORDER_COLOR));
				barFormatter = new BarFormatter(
						HIGHLIGHT_COLOR,WorkflowStatsDetailFragment.BORDER_COLOR);
			}
			else {
				barFormatter = new BarFormatter(
						fillColor,WorkflowStatsDetailFragment.BORDER_COLOR);
			}
			plot.addSeries(series, barFormatter);
			
			++workflowIndex;
		}
		List<XYSeriesRenderer> rendererList = plot.getRendererList();
		for (int i =0;i<rendererList.size();++i){
			BarRenderer tmp = (BarRenderer) rendererList.get(i);
			tmp.setBarWidth(BAR_WIDTH);
		}
		plot.setDomainBoundaries(0, workflowStatsHash.size()+1
				, BoundaryMode.FIXED);
		plot.setDomainStepValue(1);
		plot.getLegendWidget().setVisible(false);
		plot.setRangeValueFormat(NumberFormat.getIntegerInstance());
		plot.getGraphWidget().setMarginTop(GRAPH_TOP_MARGIN);
		plot.setRangeLowerBoundary(0, BoundaryMode.FIXED);
		if (maxRange == 0){
			plot.setRangeUpperBoundary(workflowIndex+1, BoundaryMode.FIXED);
		} 
		else {
			plot.setRangeUpperBoundary(maxRange 
					+ (int) plot.getRangeStepValue(), BoundaryMode.FIXED);
		}
	}
	
	/*private String shortenWorkflowName(String workflowName){
		String output = "";
		for (int i = 0; i < workflowName.length();++i){
			if (Character.isUpperCase(workflowName.charAt(i))){
				if (i==workflowName.length()-1 || Character.isUpperCase(workflowName.charAt(i+1))){
					output+=workflowName.charAt(i);
				}
				else {
					output+=workflowName.charAt(i)+".";
				}	
			}
			else if (!Character.isLetter(workflowName.charAt(i)) 
					&& !Character.isLetter(output.charAt(output.length()-1))){
						output = output.substring(0, output.length()-1) + workflowName.charAt(i);
			}
		}
		return output;
	}*/
}
