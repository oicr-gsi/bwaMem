package ca.on.oicr.pde.seqprodreporter;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import android.content.Context;
import android.view.Gravity;
import android.view.LayoutInflater;
import android.view.View;
import android.view.ViewGroup;
import android.widget.ArrayAdapter;
import android.widget.ProgressBar;
import android.widget.TextView;

public class ReportAdapter extends ArrayAdapter<Report> {

	private Context mContext;
	private List<Report> list = new ArrayList<Report>();
	private static LayoutInflater inflater = null;
	
	private static final String ESTIMATED_TIME_REMAINING_FILLER = "HH:MM";

	public ReportAdapter(Context context, int resource) {
		super(context, resource);
		this.mContext = context;
		inflater = LayoutInflater.from(mContext);
	}

	@Override
	public int getCount() {
		return this.list.size();
	}

	@Override
	public Report getItem(int position) {
		return this.list.get(position);
	}

	@Override
	public long getItemId(int position) {
		return position;
	}

	@Override
	public View getView(int position, View convertView, ViewGroup parent) {

		View newView = convertView;
		ReportHolder holder;
		Report curr = list.get(position);

		if (null == convertView) {
			holder = new ReportHolder();
			newView = inflater.inflate(R.layout.report, null);
			holder.samplename = (TextView) newView
					.findViewById(R.id.report_samplename);
			holder.wfname = (TextView) newView.findViewById(R.id.report_wfname);
			holder.ctime = (TextView) newView.findViewById(R.id.report_ctime);
			holder.lmtime = (TextView) newView.findViewById(R.id.report_lmtime);
			holder.pbar = (ProgressBar) newView.findViewById(R.id.report_pbar);
			holder.estimatedTimeRemaining = (TextView) newView.findViewById(R.id.estimated_time_remaining);
			newView.setTag(holder);
		} 
		else {
			holder = (ReportHolder) newView.getTag();
		}
		// For the case where there are no reports in the corresponding list
		// an empty report is created with a message indicating this; the text is set here
		if (list.size() == 1 && curr.getrWorkflowName().equals(Report.EMPTY_REPORT)) {
			holder.samplename.setText(curr.getrSampleName());
			holder.samplename.setGravity(Gravity.CENTER);
		}
		
		else {
			holder.samplename.setText("Sample: " + curr.getrSampleName());
			holder.wfname.setText(curr.getrWorkflowName() + " "
					+ curr.getrWorkflowVersion());
			holder.ctime.setText("Created: " + curr.getrCreateTime());
			holder.lmtime.setText("Modified: " + curr.getrLastmodTime());
			if (null != curr.getrProgress()) {
				holder.pbar.setProgress(curr.progressValue());
				holder.pbar.setVisibility(ProgressBar.VISIBLE);
				//TODO: Set the correct ETR here once back end functionality is implemented
				// Also change the text view visibility to visible in XML
				holder.estimatedTimeRemaining.setText("ETR: " + ESTIMATED_TIME_REMAINING_FILLER );
			}
			
			if (curr.getrUpdated()){
				newView.setBackgroundColor
					(getContext().getResources().getColor(R.color.highlight_green));
			}	
			else {
				newView.setBackgroundColor
					(getContext().getResources().getColor(R.color.white));
			} 
		}
		

		return newView;
	}

	static class ReportHolder {
		TextView samplename;
		TextView wfname;
		TextView ctime;
		TextView lmtime;
		ProgressBar pbar;
		TextView estimatedTimeRemaining;
	}

	public void add(Report listItem) {
		list.add(listItem);
	}

	public void sortList(Comparator<Report> comp) {
		Collections.sort(list, comp);
    }

	public void removeAllViews() {
		list.clear();
	}

}
