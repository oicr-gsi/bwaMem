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
public class WorkflowStatsDetailFragment extends Fragment {
	LinkedHashMap<String, Number[]> workflowStatsHash = new LinkedHashMap<String, Number[]>();
	String [] workflowList;
	private XYPlot completedPlot;
	private XYPlot pendingPlot;
	private XYPlot failedPlot;
	String selectedWorkflow;
	
	private boolean isTwoPane;
	
	//TODO: May need to change offset values upon further testing with different data
	private static final double TWOPANE_WORKFLOW_NAME_OFFSET = 0.18;
	private static final double TWOPANE_WORKFLOW_NUMBER_OFFSET = 0.08;
	private static final double WORKFLOW_NAME_OFFSET = 0.25;
	private static final double WORKFLOW_NUMBER_OFFSET = 0.09;
	
	private static final float GRAPH_TOP_MARGIN = 15f;
	private static final float BAR_WIDTH = 40f;
	
	private static final int BORDER_COLOR = 0xFF000000;
	private static final int TRANSPARENT_COLOR = 0x30000000;
	private static final int HIGHLIGHT_COLOR = 0xFFFFFFFF;
	private static final int FILL_COLOR_COMPLETED = 0XFF00FF00 ;
	private static final int FILL_COLOR_PENDING = 0xFFFFFF00;
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
		isTwoPane = false;
		if (getArguments().containsKey(ARG_ITEM_ID) && getArguments().containsKey("WorkflowList")) {
			selectedWorkflow = getArguments().getString(ARG_ITEM_ID);
			workflowList = getArguments().getStringArray("WorkflowList");
			for (int i = 0;i<workflowList.length; ++i){
				getWorkflowStats(workflowList[i]);
			}
			if (getArguments().containsKey("IsTwoPane")){
				isTwoPane = true;
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
		
		setUpXYPlot(completedPlot, 0);
		setUpXYPlot(pendingPlot, 1);
		setUpXYPlot(failedPlot, 2);
		
		return rootView;
	}
	
	private void setUpXYPlot(XYPlot plot, int workflowType){
		int workflowIndex = 1;
		int fillColor;
		int maxRange = 0;
	
		switch(workflowType){
			case 0:
				fillColor = FILL_COLOR_COMPLETED;
				break;
			case 1:
				fillColor = FILL_COLOR_FAILED;
				break;
			case 2: 
				fillColor = FILL_COLOR_PENDING;
				break;
			default:
				fillColor = BORDER_COLOR;
				break;
		}
				
		for (String workflowName : workflowStatsHash.keySet()){
			
			if (1 == workflowIndex){
				maxRange = workflowStatsHash.get(workflowName)[workflowType].intValue();
			}
			else if (maxRange < workflowStatsHash.get(workflowName)[workflowType].intValue()){
				maxRange = workflowStatsHash.get(workflowName)[workflowType].intValue();
			}
			
			BarFormatter barFormatter;
			SimpleXYSeries series = new SimpleXYSeries(workflowName);
			series.addFirst(workflowIndex, 
					workflowStatsHash.get(workflowName)[workflowType]);
			
			String workflowNameEdit = workflowName.length() < 11 ? workflowName : shortenWorkflowName(workflowName);
			double nameOffset = isTwoPane ? TWOPANE_WORKFLOW_NAME_OFFSET : WORKFLOW_NAME_OFFSET;
			double numberOffset = isTwoPane ? TWOPANE_WORKFLOW_NUMBER_OFFSET : WORKFLOW_NUMBER_OFFSET;
			
			if (!selectedWorkflow.equals("NoSelectedWorkflow") 
					&& workflowName.equals(selectedWorkflow)){
				plot.addMarker(new XValueMarker(workflowIndex - nameOffset,workflowNameEdit,
						new YPositionMetric(0,YLayoutStyle.ABSOLUTE_FROM_TOP),
						0,TRANSPARENT_COLOR));
				plot.addMarker(new XValueMarker(workflowIndex- numberOffset,
						workflowStatsHash.get(workflowName)[workflowType].toString(),
							new YPositionMetric(0,YLayoutStyle.ABSOLUTE_FROM_BOTTOM),
								0,BORDER_COLOR));
				barFormatter = new BarFormatter(
						HIGHLIGHT_COLOR,BORDER_COLOR);
			}
			else {
				plot.addMarker(new XValueMarker(workflowIndex - nameOffset,workflowNameEdit,
						new YPositionMetric(0,YLayoutStyle.ABSOLUTE_FROM_TOP),
							0,BORDER_COLOR));
				barFormatter = new BarFormatter(
						fillColor,BORDER_COLOR);
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
		if (maxRange == 0){
			plot.setRangeUpperBoundary(workflowIndex+1, BoundaryMode.FIXED);
		} 
		else {
			plot.setRangeUpperBoundary(maxRange 
					+ (int) plot.getRangeStepValue(), BoundaryMode.FIXED);
		}
	}
	private String shortenWorkflowName(String workflowName){
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
	}
}
