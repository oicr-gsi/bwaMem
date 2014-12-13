package ca.on.oicr.pde.seqprodreporter;

import android.app.Activity;
import android.app.FragmentManager;
import android.app.FragmentTransaction;
import android.os.Bundle;
import android.util.Log;



public class WorkflowStatsActivity extends Activity {
    
	private FragmentManager mFragmentManager;
	
	protected final static  String TAG = "Reporter Stats";
	
	@Override
	protected void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		
		Log.i(TAG, getClass().getSimpleName() + ":entered onCreate()");
		setContentView(R.layout.activity_workflow_stats);

		this.mFragmentManager = getFragmentManager();
		FragmentTransaction fragmentTransaction = mFragmentManager
				.beginTransaction();
		if (this.isLayoutLarge()) {
			fragmentTransaction.add(R.id.piechart_container,
					new WorkflowChartContainerFragment());
			fragmentTransaction.add(R.id.workflow_list_container,
					new WorkflowListFragment());
		} else {	
			fragmentTransaction.add(R.id.fragment_pager,
				new WorkflowListFragment());
		}
		fragmentTransaction.commit();

	}

	private boolean isLayoutLarge() {
		return findViewById(R.id.fragment_pager) == null;
	}
	
}
