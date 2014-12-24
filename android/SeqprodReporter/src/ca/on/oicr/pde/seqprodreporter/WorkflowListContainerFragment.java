package ca.on.oicr.pde.seqprodreporter;

import android.app.Fragment;
import android.app.FragmentTransaction;
import android.content.Intent;
import android.os.Bundle;
import android.view.LayoutInflater;
import android.view.View;
import android.view.View.OnClickListener;
import android.view.ViewGroup;
import android.widget.TextView;

/**
 * A list fragment representing a list of WorkflowStats. This fragment also
 * supports tablet devices by allowing list items to be given an 'activated'
 * state upon selection. This helps indicate which item is currently being
 * viewed in a {@link WorkflowStatsDetailFragment}.
 * <p>
 * Activities containing this fragment MUST implement the {@link OnItemSelectedListener}
 * interface.
 */

/**
 * @author pruzanov
 * 
 * Re-structuring modules for more rational code organisation
 */

public class WorkflowListContainerFragment extends Fragment {
	
	private TextView completedTextView;
	private TextView failedTextView;
	private TextView pendingTextView;
	//private WorkflowStatsListFragment listFragment;
    private int[] workflowTypesTotal;
    private String[] activeWorkflows;
	
    private String[] getActiveWorkflows() {
		return activeWorkflows;
	}

	private void setActiveWorkflows(String[] activeWorkflows) {
		this.activeWorkflows = activeWorkflows;
	}

	private int getWorkflowTypesTotal(int index) {
    	if (index >= 0 && index <= this.workflowTypesTotal.length - 1)
		    return workflowTypesTotal[index];
    	else
    		return 0;
	}

	protected void setWorkflowTypesTotal(int[] workflowTypesTotal) {
		this.workflowTypesTotal = workflowTypesTotal;
	}

	/**
	 * Mandatory empty constructor for the fragment manager to instantiate the
	 * fragment (e.g. upon screen orientation changes).
	 */
	public WorkflowListContainerFragment() {
	}
	
	private OnClickListener onTextViewClick = new OnClickListener() {
		@Override
		public void onClick(View textView) {
			Intent intent = new Intent(getActivity(),
					ReporterActivity.class);
			intent.setAction(Intent.ACTION_RUN);
			if (textView.equals(completedTextView)) {
				intent.putExtra("selectedTab",
						ReporterActivity.COMPLETED_WORKFLOW_TAB_INDEX);
			} else if (textView.equals(failedTextView)) {
				intent.putExtra("selectedTab",
						ReporterActivity.FAILED_WORKFLOW_TAB_INDEX);
			} else if (textView.equals(pendingTextView)) {
				intent.putExtra("selectedTab",
						ReporterActivity.PENDING_WORKFLOW_TAB_INDEX);
			}
			startActivity(intent);
		}
	};

	public static WorkflowListContainerFragment InstanceOf(int[] typesTotal, String[] workflows) {
		WorkflowListContainerFragment fragment = new WorkflowListContainerFragment();
		fragment.setWorkflowTypesTotal(typesTotal);
		fragment.setActiveWorkflows(workflows);
		return fragment;
	}
	
	@Override 
	public void onCreate(Bundle savedInstanceState){
		super.onCreate(savedInstanceState);
	}

	@Override
	public View onCreateView(LayoutInflater inflater, ViewGroup container,
			Bundle savedInstanceState) {
	
		View rootView = inflater.inflate(
				R.layout.fragment_workflowstats_container, container, false);
		
		//Setup Text fields
		this.completedTextView = (TextView) rootView.findViewById(R.id.total_number_of_completed);
		completedTextView.setText("Completed: " + getWorkflowTypesTotal(0));
		completedTextView.setOnClickListener(onTextViewClick);

		this.failedTextView = (TextView) rootView.findViewById(R.id.total_number_of_failed);
		failedTextView.setText("Failed: " + getWorkflowTypesTotal(1));
		failedTextView.setOnClickListener(onTextViewClick);

		this.pendingTextView = (TextView) rootView.findViewById(R.id.total_number_of_pending);
		pendingTextView.setText("Pending: " + getWorkflowTypesTotal(2));
		pendingTextView.setOnClickListener(onTextViewClick);
		
		//TODO setup list fragment
		FragmentTransaction fragmentTransaction = getFragmentManager()
				.beginTransaction();

		fragmentTransaction.add(R.id.workflow_list_container,
					            WorkflowListFragment.InstanceOf(this.getActiveWorkflows()));
		fragmentTransaction.commit();
		return rootView;
	}


}
