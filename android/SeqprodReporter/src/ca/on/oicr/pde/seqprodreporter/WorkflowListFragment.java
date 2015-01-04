package ca.on.oicr.pde.seqprodreporter;

import android.app.Activity;
import android.app.ListFragment;
import android.content.Intent;
import android.os.Bundle;
import android.util.Log;
import android.view.LayoutInflater;
import android.view.View;
import android.view.View.OnClickListener;
import android.view.ViewGroup;
import android.widget.ArrayAdapter;
import android.widget.ListView;
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
 *         Re-structuring modules for more rational code organisation
 */

public class WorkflowListFragment extends ListFragment {

	private TextView completedTextView;
	private TextView failedTextView;
	private TextView pendingTextView;

	private int[] workflowTypesTotal;
	private String[] workflowNames;
	private ArrayAdapter<String> mAdapter;
	/**
	 * The serialization (saved instance state) Bundle key representing the
	 * activated item position. Only used on tablets.
	 */
    public static final String ALL_WORKFLOWS = "All Workflows";
	/**
	 * The fragment's current callback object, which is notified of list item
	 * clicks.
	 */
	private OnItemSelectedListener mCallbacks = sDummyCallbacks;

	/**
	 * A callback interface that all activities containing this fragment must
	 * implement. This mechanism allows activities to be notified of item
	 * selections.
	 */
	public interface OnItemSelectedListener {
		/**
		 * Callback for when an item has been selected.
		 */
		public void onItemSelected(String id);
	}

	private OnClickListener onTextViewClick = new OnClickListener() {
		@Override
		public void onClick(View textView) {
			Intent intent = new Intent(getActivity(), ReporterActivity.class);
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
			Log.v(ReporterActivity.TAG, "Clicked on text "
					+ ((TextView) textView).getText());
		}
	};

	/**
	 * A dummy implementation of the {@link OnItemSelectedListener} interface
	 * that does nothing. Used only when this fragment is not attached to an
	 * activity.
	 */
	private static OnItemSelectedListener sDummyCallbacks = new OnItemSelectedListener() {
		@Override
		public void onItemSelected(String id) {
		}
	};

	private int getWorkflowTypesTotal(int index) {
		if (index >= 0 && index <= this.workflowTypesTotal.length - 1)
			return workflowTypesTotal[index];
		else
			return 0;
	}

	protected void setWorkflowTypesTotal(int[] workflowTypesTotal) {
		this.workflowTypesTotal = workflowTypesTotal;
	}

	private void setWorkflowNames(String[] workflowNames) {
		this.workflowNames = workflowNames;
	}
	/**
	 * Mandatory empty constructor for the fragment manager to instantiate the
	 * fragment (e.g. upon screen orientation changes).
	 */
	public WorkflowListFragment() {
	}

	// TODO add grand totals
	public static WorkflowListFragment InstanceOf(int[] totals,
			String[] workflows) {
		WorkflowListFragment fragment = new WorkflowListFragment();
		fragment.setWorkflowNames(workflows);
		fragment.setWorkflowTypesTotal(totals);
		return fragment;
	}

	@Override
	public void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		mAdapter = new ArrayAdapter<String>(getActivity(),
				android.R.layout.simple_list_item_activated_1,
				this.workflowNames);
		setListAdapter(mAdapter);
	}

	@Override
	public View onCreateView(LayoutInflater inflater, ViewGroup container,
			Bundle savedInstanceState) {

		View rootView = inflater.inflate(
				R.layout.fragment_workflowstats_container, container, false);

		// Setup Text fields
		this.completedTextView = (TextView) rootView.findViewById(R.id.total_number_of_completed);
		completedTextView.setText("Completed: " + getWorkflowTypesTotal(0));
		completedTextView.setOnClickListener(onTextViewClick);

		this.failedTextView = (TextView) rootView.findViewById(R.id.total_number_of_failed);
		failedTextView.setText("Failed: " + getWorkflowTypesTotal(1));
		failedTextView.setOnClickListener(onTextViewClick);

		this.pendingTextView = (TextView) rootView.findViewById(R.id.total_number_of_pending);
		pendingTextView.setText("Pending: " + getWorkflowTypesTotal(2));
		pendingTextView.setOnClickListener(onTextViewClick);

		return rootView;
	}

	@Override
	public void onAttach(Activity activity) {
		super.onAttach(activity);
		// Activities containing this fragment must implement its callbacks.
		if (!(activity instanceof OnItemSelectedListener)) {
			throw new IllegalStateException(
					"Activity must implement fragment's callbacks.");
		}

		mCallbacks = (OnItemSelectedListener) activity;
	}

	@Override
	public void onDetach() {
		super.onDetach();
		// Reset the active callbacks interface to the dummy implementation.
		mCallbacks = sDummyCallbacks;
	}

	@Override
	public void onListItemClick(ListView listView, View view, int position,
			long id) {

		String workflowName = (String) ((TextView) view).getText();
		mCallbacks.onItemSelected(workflowName);
	}

}
