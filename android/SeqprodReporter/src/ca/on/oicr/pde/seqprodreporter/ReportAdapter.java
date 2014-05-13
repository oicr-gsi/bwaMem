package ca.on.oicr.pde.seqprodreporter;

import java.util.ArrayList;

import android.R.color;
import android.content.Context;
import android.view.LayoutInflater;
import android.view.View;
import android.view.ViewGroup;
import android.widget.BaseAdapter;
import android.widget.ProgressBar;
import android.widget.TextView;

public class ReportAdapter extends BaseAdapter {
	private Context mContext;
	//TODO (PDE-577) need to flag updated reports, for this we need to extract time and compare it to the time of the last update
	//This parameter may be local to this adapter
	private ArrayList<Report> list = new ArrayList<Report>();
	private static LayoutInflater inflater = null;
	private final int UPDATE_COLOR = color.holo_orange_light;
		
	public ReportAdapter(Context c) {
		this.mContext = c;
		inflater = LayoutInflater.from(mContext);
	}

	@Override
	public int getCount() {
		return this.list.size();
	}

	@Override
	public Object getItem(int position) {
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
			holder.samplename = (TextView) newView.findViewById(R.id.report_samplename);
			holder.wfname = (TextView) newView.findViewById(R.id.report_wfname);
			holder.ctime  = (TextView) newView.findViewById(R.id.report_ctime);
			holder.lmtime = (TextView) newView.findViewById(R.id.report_lmtime);
			holder.pbar   = (ProgressBar) newView.findViewById(R.id.report_pbar);
			newView.setTag(holder);
		} else {
			holder = (ReportHolder) newView.getTag();
		}

		holder.samplename.setText("Sample: " + curr.getrSampleName());
		holder.wfname.setText(curr.getrWorkflowName() + " " + curr.getrWorkflowVersion());
		holder.ctime.setText("Created: " + curr.getrCreateTime());
		holder.lmtime.setText("Modified: " + curr.getrLastmodTime());
		if (null != curr.getrProgress()) {
			holder.pbar.setProgress(curr.progressValue());
			holder.pbar.setVisibility(ProgressBar.VISIBLE);	
		}
		if (curr.getrUpdated())
			//TODO (PDE-577) this doesn't have any effect at the moment, need to debug
			// perhaps, invalidate the corresponding View?
			holder.samplename.setBackgroundColor(UPDATE_COLOR);

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
		notifyDataSetChanged();
	}
	
	public ArrayList<Report> getList(){
		return list;
	}
	
	public void removeAllViews(){
		list.clear();
		this.notifyDataSetChanged();
	}

}
