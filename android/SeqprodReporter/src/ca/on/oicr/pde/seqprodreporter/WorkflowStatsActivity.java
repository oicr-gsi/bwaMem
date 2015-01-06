package ca.on.oicr.pde.seqprodreporter;

import android.app.Activity;
import android.app.FragmentManager;
import android.app.FragmentTransaction;
import android.content.SharedPreferences;
import android.database.Cursor;
import android.os.Bundle;
import android.util.Log;
import ca.on.oicr.pde.seqprodprovider.DataContract;

/**
 * @author pruzanov@oicr.on.ca
 * 
 * Re-structuring modules for more rational code organization
 */

public class WorkflowStatsActivity extends Activity implements
		WorkflowListFragment.OnItemSelectedListener {

	private FragmentManager mFragmentManager;
	private WorkflowChartFragment pieChartFragment;
	protected final static String ALL_WORKFLOWS = "All Workflows";
	protected final static String WORKFLOW_PIE_CHART_VALUES = "WorkflowTypeTotals";
	protected final static String WORKFLOW_SELECTED = "SelectedWorkflow";
	protected final static String CHART_SHOWN = "ChartViewActive";
	protected final static String TAG = "Reporter Stats";

	protected String[] activeWorkflows;
	protected int[] selectedTotals;
	protected int[] grandTotals;
	private String selectedWorkflow;
	private boolean chartViewActive;

	@Override
	protected void onCreate(Bundle savedInstanceState) {

		super.onCreate(savedInstanceState);
		Log.i(TAG, getClass().getSimpleName() + ":entered onCreate()");
		setContentView(R.layout.activity_workflow_stats);

		this.activeWorkflows = this.getWorkflows();
		this.selectedTotals  = this.getPieChartValues(this.selectedWorkflow);
		this.chartViewActive = false;

		this.mFragmentManager = getFragmentManager();
		FragmentTransaction fragmentTransaction = mFragmentManager
				.beginTransaction();
		WorkflowListFragment listFragment = WorkflowListFragment.InstanceOf(
				this.selectedTotals, this.activeWorkflows);

		if (this.isLayoutLarge()) {
			this.pieChartFragment = WorkflowChartFragment.InstanceOf(
					this.selectedTotals, this.selectedWorkflow);
			fragmentTransaction.add(R.id.piechart_container,
					this.pieChartFragment);
			fragmentTransaction.add(R.id.workflow_list_container, listFragment);
		} else {
			fragmentTransaction.add(R.id.fragment_pager, listFragment);
		}
		fragmentTransaction.commit();
	}

	private boolean isLayoutLarge() {
		return findViewById(R.id.fragment_pager) == null;
	}

	@Override
	public void onRestoreInstanceState(Bundle savedInstanceState) {
		if (null != savedInstanceState) {
			try {
				this.selectedTotals = savedInstanceState
						.getIntArray(WORKFLOW_PIE_CHART_VALUES);
				this.selectedWorkflow = savedInstanceState
						.getString(WORKFLOW_SELECTED);
				this.chartViewActive = savedInstanceState.getBoolean(
						CHART_SHOWN, false);
				this.refreshChartView();
			} catch (Exception e) {
				Log.d(ReporterActivity.TAG,
						"Workflow totals could not be retrieved");
			}
		}
		super.onRestoreInstanceState(savedInstanceState);
	}

	@Override
	public void onSaveInstanceState(Bundle outState) {
		if (null != this.selectedTotals)
			outState.putIntArray(WORKFLOW_PIE_CHART_VALUES, this.selectedTotals);
		outState.putString(WORKFLOW_SELECTED, this.selectedWorkflow);
		outState.putBoolean(CHART_SHOWN, this.chartViewActive);
	}

	@Override
	public void onItemSelected(String id) {

		Log.d(TAG, "Received Click from " + id);
		this.selectedWorkflow = id;
		this.selectedTotals = this.getPieChartValues(id);
		this.chartViewActive = true;
		this.refreshChartView();

	}

	private void refreshChartView() {

		if (this.isLayoutLarge()) {
			// Just update pie chart
			this.pieChartFragment.updatePieChartValues(this.selectedTotals,
					                                   this.selectedWorkflow);
		} else {
			// We received a click, which means pie chart is not visible. Need to
			// create new pie chart fragment and replace the currently shown one
			if (this.chartViewActive) {
				FragmentTransaction pieChartShowTransaction = mFragmentManager
						.beginTransaction()
						.addToBackStack(CHART_SHOWN)
						.replace(R.id.fragment_pager,
								 WorkflowChartFragment.InstanceOf(this.selectedTotals,
										                          this.selectedWorkflow),
						         "Chart");
				pieChartShowTransaction.commit();
			}
		}
	}

	/**
	 * function for updating pie chart (all MySQL code needs to be here)
	 * 
	 * @return workflow names as String array
	 */
	private String[] getWorkflows() {
		SharedPreferences sp = getSharedPreferences(
				ReporterActivity.PREFERENCE_FILE, ReporterActivity.MODE_PRIVATE);
		String timeRange = sp.getString("pref_summaryScope", getResources()
				.getStringArray(R.array.pref_summaryScope_entries)[0]);

		long earliest = ReporterActivity.getEarliestMillis(timeRange);
		Cursor c = this.getContentResolver().query(DataContract.CONTENT_URI,
				new String[] { "DISTINCT " + DataContract.WORKFLOW },
				DataContract.LM_TIME + "> ? ", new String[] { earliest + "" },
				DataContract.WORKFLOW + " ASC");

		String[] workflowNames = new String[c.getCount() + 1];
		workflowNames[0] = ALL_WORKFLOWS;
		int index = 1;
		if (c.moveToFirst()) {
			do {
				workflowNames[index] = c.getString(
						c.getColumnIndex(DataContract.WORKFLOW));
				++index;
			} while (c.moveToNext());
		}
		c.close();
		return workflowNames;
	}

	/**
	 * function for updating totals to be used with ChartWidget can have null as
	 * argument to get totals for list of all wfs
	 * 
	 * @param selectedWorkflow
	 * @return
	 */
	private int[] getPieChartValues(String selectedWorkflow) {
		int[] selectedWorkflowNumbers = new int[ReporterActivity.types.length];
		SharedPreferences sp = getSharedPreferences(
				ReporterActivity.PREFERENCE_FILE, ReporterActivity.MODE_PRIVATE);
		String timeRange = sp.getString("pref_summaryScope", getResources()
				.getStringArray(R.array.pref_summaryScope_entries)[0]);

		long earliest = ReporterActivity.getEarliestMillis(timeRange);
		for (int i = 0; i < ReporterActivity.types.length; ++i) {
			String[] projection;
			StringBuilder selection = new StringBuilder();
			String[] selectionArgs;
			if (selectedWorkflow == null
					|| selectedWorkflow.equals(ALL_WORKFLOWS)) {
				projection = new String[] { DataContract.WR_TYPE };
				selectionArgs = new String[] { ReporterActivity.types[i],
						"" + earliest };
			} else {
				projection = new String[] { DataContract.WORKFLOW };
				selectionArgs = new String[] { selectedWorkflow,
						ReporterActivity.types[i], "" + earliest };
				selection.append(DataContract.WORKFLOW + "=? AND ");
			}
			selection.append(DataContract.WR_TYPE + "=? AND "
					+ DataContract.LM_TIME + "> ? ");

			Cursor c = this.getContentResolver().query(
					DataContract.CONTENT_URI, projection, selection.toString(),
					selectionArgs, null);

			selectedWorkflowNumbers[i] = c.getCount();
			c.close();
		}

		return selectedWorkflowNumbers;
	}

}
