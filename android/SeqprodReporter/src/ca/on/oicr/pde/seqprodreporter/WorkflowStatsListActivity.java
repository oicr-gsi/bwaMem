package ca.on.oicr.pde.seqprodreporter;

import java.util.ArrayList;

import ca.on.oicr.pde.seqprodprovider.DataContract;
import android.content.Intent;
import android.database.Cursor;
import android.os.Bundle;
import android.support.v4.app.FragmentActivity;
import android.support.v4.app.NavUtils;
import android.view.MenuItem;
import android.view.View;
import android.view.View.OnClickListener;
import android.widget.TextView;

/**
 * An activity representing a list of WorkflowStats. This activity has different
 * presentations for handset and tablet-size devices. On handsets, the activity
 * presents a list of items, which when touched, lead to a
 * {@link WorkflowStatsDetailActivity} representing item details. On tablets,
 * the activity presents the list of items and item details side-by-side using
 * two vertical panes.
 * <p>
 * The activity makes heavy use of fragments. The list of items is a
 * {@link WorkflowStatsListFragment} and the item details (if present) is a
 * {@link WorkflowStatsDetailFragment}.
 * <p>
 * This activity also implements the required
 * {@link WorkflowStatsListFragment.Callbacks} interface to listen for item
 * selections.
 */
public class WorkflowStatsListActivity extends FragmentActivity implements
		WorkflowStatsListFragment.Callbacks {
	private TextView completedTextView;
	private TextView failedTextView;
	private TextView pendingTextView;
	
	private WorkflowStatsListFragment listFragment;
	
	/**
	 * Whether or not the activity is in two-pane mode, i.e. running on a tablet
	 * device.
	 */
	private boolean mTwoPane;

	@Override
	protected void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		setContentView(R.layout.activity_workflowstats_list);
		// Show the Up button in the action bar.
		getActionBar().setDisplayHomeAsUpEnabled(true);
		
		listFragment = (WorkflowStatsListFragment) getSupportFragmentManager().findFragmentById(R.id.workflowstats_list);
		setUpTextViews();
		
		if (findViewById(R.id.workflowstats_detail_container) != null) {
			// The detail container view will be present only in the
			// large-screen layouts (res/values-large and
			// res/values-sw600dp). If this view is present, then the
			// activity should be in two-pane mode.
			mTwoPane = true;
			replaceDetailsFragment("NoSelectedWorkflow");
			// In two-pane mode, list items should be given the
			// 'activated' state when touched.
			((WorkflowStatsListFragment) getSupportFragmentManager()
					.findFragmentById(R.id.workflowstats_list))
					.setActivateOnItemClick(true);
		}

		// TODO: If exposing deep links into your app, handle intents here.
	}
	
	private OnClickListener onTextViewClick = new OnClickListener(){
		@Override
		public void onClick(View textView){
			Intent intent = new Intent(WorkflowStatsListActivity.this, ReporterActivity.class);
			if (textView.equals(completedTextView)){
				intent.putExtra("selectedTab", ReporterActivity.COMPLETED_WORKFLOW_TAB_INDEX);
			}
			else if (textView.equals(failedTextView)){
				intent.putExtra("selectedTab", ReporterActivity.FAILED_WORKFLOW_TAB_INDEX);
			}
			else if (textView.equals(pendingTextView)){
				intent.putExtra("selectedTab", ReporterActivity.PENDING_WORKFLOW_TAB_INDEX);
			}
			startActivity(intent);
		} 
		
		
	};
	
	private int workflowTypeTotal(String workflowType){
		Cursor c = this.getContentResolver()
				.query(DataContract.CONTENT_URI, new String[]{DataContract.WR_TYPE}, DataContract.WR_TYPE + "=?", new String[]{workflowType} , null);			
		return c.getCount();
	}
	
	private void setUpTextViews(){
		this.completedTextView = (TextView) findViewById(R.id.total_number_of_completed);
		completedTextView.setText("Completed: " + workflowTypeTotal(ReporterActivity.types[0]));
		completedTextView.setOnClickListener(onTextViewClick);
		
		this.failedTextView = (TextView) findViewById(R.id.total_number_of_failed);
		failedTextView.setText("Failed: " + workflowTypeTotal(ReporterActivity.types[1]));
		failedTextView.setOnClickListener(onTextViewClick);
		
		this.pendingTextView = (TextView) findViewById(R.id.total_number_of_pending);
		pendingTextView.setText("Pending: " + workflowTypeTotal(ReporterActivity.types[2]));
		pendingTextView.setOnClickListener(onTextViewClick);
	}

	@Override
	public boolean onOptionsItemSelected(MenuItem item) {
		int id = item.getItemId();
		if (id == android.R.id.home) {
			// This ID represents the Home or Up button. In the case of this
			// activity, the Up button is shown. Use NavUtils to allow users
			// to navigate up one level in the application structure. For
			// more details, see the Navigation pattern on Android Design:
			//
			// http://developer.android.com/design/patterns/navigation.html#up-vs-back
			//
			NavUtils.navigateUpFromSameTask(this);
			return true;
		}
		return super.onOptionsItemSelected(item);
	}

	/**
	 * Callback method from {@link WorkflowStatsListFragment.Callbacks}
	 * indicating that the item with the given ID was selected.
	 */
	@Override
	public void onItemSelected(String id) {

		if (mTwoPane) {
			// In two-pane mode, show the detail view in this activity by
			// adding or replacing the detail fragment using a
			// fragment transaction.
			replaceDetailsFragment(id);

		} else {
			// In single-pane mode, simply start the detail activity
			// for the selected item ID.
			Intent detailIntent = new Intent(this,
					WorkflowStatsDetailActivity.class);
			detailIntent.putExtra("WorkflowList",listFragment.getWorkflowList()
					.toArray(new String[listFragment.getWorkflowList().size()]));
			detailIntent.putExtra(WorkflowStatsDetailFragment.ARG_ITEM_ID, id);
			startActivity(detailIntent);
		}
	}
	
	private void replaceDetailsFragment(String id){
		Bundle arguments = new Bundle();
		arguments.putStringArray("WorkflowList", listFragment.getWorkflowList()
				.toArray(new String[listFragment.getWorkflowList().size()]));
		arguments.putString(WorkflowStatsDetailFragment.ARG_ITEM_ID, id);
		WorkflowStatsDetailFragment fragment = new WorkflowStatsDetailFragment();
		fragment.setArguments(arguments);
		getSupportFragmentManager().beginTransaction()
				.replace(R.id.workflowstats_detail_container, fragment)
				.commit();
	}
}
