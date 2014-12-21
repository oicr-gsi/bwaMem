package ca.on.oicr.pde.seqprodreporter;

import android.app.Activity;
import android.app.FragmentManager;
import android.app.FragmentTransaction;
import android.content.SharedPreferences;
import android.database.Cursor;
import android.os.Bundle;
import android.util.Log;
import ca.on.oicr.pde.seqprodprovider.DataContract;



public class WorkflowStatsActivity extends Activity implements
    WorkflowListFragment.OnItemSelectedListener {
    
	private FragmentManager mFragmentManager;
	protected final static String ALL_WORKFLOWS = "All Workflows";
	protected final static String TAG = "Reporter Stats";
	
	@Override
	protected void onCreate(Bundle savedInstanceState) {
		// TODO get values from db
		super.onCreate(savedInstanceState);
		
		Log.i(TAG, getClass().getSimpleName() + ":entered onCreate()");
		setContentView(R.layout.activity_workflow_stats);

		this.mFragmentManager = getFragmentManager();
		FragmentTransaction fragmentTransaction = mFragmentManager
				.beginTransaction();
		if (this.isLayoutLarge()) {
			fragmentTransaction.add(R.id.piechart_container,
					new WorkflowChartFragment());
			fragmentTransaction.add(R.id.workflow_list_container,
					new WorkflowListFragment());
		} else {	
			fragmentTransaction.add(R.id.fragment_pager,
				new WorkflowListFragment());
		}
		fragmentTransaction.commit();
		//TODO update values for PieChart and List

	}

	private boolean isLayoutLarge() {
		return findViewById(R.id.fragment_pager) == null;
	}

	@Override
	public void onItemSelected(String id) {
		// TODO Auto-generated method stub
		// React on list item selected
		
	}
	
	//function for updating piechart (all MySQL code needs to be here)
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
	
	//function for updating totals to be used with ChartWidget
	private int[] getPieChartValues(String selectedWorkflow){
		int [] selectedWorkflowNumbers = new int[ReporterActivity.types.length];
		SharedPreferences sp = getSharedPreferences(ReporterActivity.PREFERENCE_FILE, 
				                                    ReporterActivity.MODE_PRIVATE);	
		String timeRange = sp.getString("pref_summaryScope", getResources().getStringArray(
				R.array.pref_summaryScope_entries)[0]);

		long earliest = ReporterActivity.getEarliestMillis(timeRange);
			for (int i = 0; i < ReporterActivity.types.length; ++i){
				Cursor c = this.getContentResolver().query(
						DataContract.CONTENT_URI,
						new String[]{DataContract.WORKFLOW},
						             DataContract.WORKFLOW + "=? AND " 
						           + DataContract.WR_TYPE + "=? AND "
				                   + DataContract.LM_TIME + "> ? ",
						new String[]{selectedWorkflow,
								     ReporterActivity.types[i],
								     "" + earliest},
						null);
				
				selectedWorkflowNumbers[i] = c.getCount();
				c.close();
			}
		
		return selectedWorkflowNumbers;
	}
}
