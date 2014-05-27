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
	// TODO (PDE-577) need to flag updated reports, for this we need to extract
	// time and compare it to the time of the last update
	// This parameter may be local to this adapter
	private List<Report> list = new ArrayList<Report>();
	private static LayoutInflater inflater = null;
	private final int UPDATE_COLOR = 0xFFCCFF99;
	private final int DEFAULT_COLOR = 0xFFFFFFFF;
	

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
			}
			//TODO: Set color values as static values
	
			if (curr.getrUpdated()){
				newView.setBackgroundColor(UPDATE_COLOR);
			}	
			else {
				newView.setBackgroundColor(DEFAULT_COLOR);
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
