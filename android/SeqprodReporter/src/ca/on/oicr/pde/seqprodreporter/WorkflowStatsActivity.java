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
 * @author pruzanov
 * 
 * Re-structuring modules for more rational code organisation
 */

public class WorkflowStatsActivity extends Activity implements
    WorkflowListFragment.OnItemSelectedListener {
    
	private FragmentManager mFragmentManager;
	protected final static String ALL_WORKFLOWS = "All Workflows";
	protected final static String TAG = "Reporter Stats";
	// List/Chart values, may be requested multiple times
	protected String[] activeWorkflows;
	protected int[] selectedTotals;
	protected int[] grandTotals;
	
	@Override
	protected void onCreate(Bundle savedInstanceState) {
		
		super.onCreate(savedInstanceState);
		Log.i(TAG, getClass().getSimpleName() + ":entered onCreate()");
		setContentView(R.layout.activity_workflow_stats);
        // TODO get values from db
		this.mFragmentManager = getFragmentManager();
		FragmentTransaction fragmentTransaction = mFragmentManager
				.beginTransaction();
		if (this.isLayoutLarge()) {
			fragmentTransaction.add(R.id.piechart_container,
					WorkflowChartFragment.InstanceOf(this.selectedTotals));
			fragmentTransaction.add(R.id.workflow_list_container,
					//TODO add grand totals
					WorkflowListFragment.InstanceOf(this.activeWorkflows));
		} else {	
			fragmentTransaction.add(R.id.fragment_pager,
				new WorkflowListFragment());
		}
		fragmentTransaction.commit();

	}

	private boolean isLayoutLarge() {
		return findViewById(R.id.fragment_pager) == null;
	}

	@Override
	public void onItemSelected(String id) {
		// TODO Auto-generated method stub
		// React on list item selected
		// chart.updateTypeData(newValues, wfName);
		// OR replace fragment and update
		
	}
	
	// function for updating pie chart (all MySQL code needs to be here)
	private String[] getWorkflows() {
		Cursor c = this.getContentResolver()
				  .query(DataContract.CONTENT_URI,
					     new String[]{"DISTINCT " + DataContract.WORKFLOW},
					     null, 
					     null,
					     DataContract.WORKFLOW + " ASC");
		
		String [] workflowNames = new String[c.getCount() + 1];
		workflowNames[0] = ALL_WORKFLOWS;
		int index = 1;
		if (c.moveToFirst()){
			do {
				workflowNames[index] = c.getString(c.getColumnIndex(DataContract.WORKFLOW));
				++index;
			} while (c.moveToNext());
		}
		c.close();
		return workflowNames;
	}
	
	// function for updating totals to be used with ChartWidget
	// can have null as argument to get totals for list of all wfs
	private int[] getPieChartValues(String selectedWorkflow){
		int [] selectedWorkflowNumbers = new int[ReporterActivity.types.length];
		SharedPreferences sp = getSharedPreferences(ReporterActivity.PREFERENCE_FILE, 
				                                    ReporterActivity.MODE_PRIVATE);	
		String timeRange = sp.getString("pref_summaryScope", getResources().getStringArray(
				R.array.pref_summaryScope_entries)[0]);

		long earliest = ReporterActivity.getEarliestMillis(timeRange);
			for (int i = 0; i < ReporterActivity.types.length; ++i){
				String[] projection;
				StringBuilder selection = new StringBuilder();
				String[] selectionArgs;
				if (selectedWorkflow == null) {
					projection = new String[] { DataContract.WR_TYPE };
					selectionArgs = new String[] {ReporterActivity.types[i], 
							                      "" + earliest};
				} else {
					projection = new String[] {DataContract.WORKFLOW};
					selectionArgs = new String[] {selectedWorkflow,
						                          ReporterActivity.types[i],
						                          "" + earliest };
					selection.append(DataContract.WORKFLOW + "=? AND ");
				}
				selection.append(DataContract.WR_TYPE + "=? AND "
				               + DataContract.LM_TIME + "> ? ");
				
				Cursor c = this.getContentResolver().query(
			               DataContract.CONTENT_URI,
			               projection,
			               selection.toString(),
			               selectionArgs,
			               null);
				
				selectedWorkflowNumbers[i] = c.getCount();
				c.close();
			}
		
		return selectedWorkflowNumbers;
	}
}
