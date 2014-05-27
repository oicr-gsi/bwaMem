package ca.on.oicr.pde.seqprodreporter;

import android.text.format.Time;
import android.util.Log;
import android.util.TimeFormatException;

public class Report {
	private String rSampleName;
	private String rWorkflowName;
	private String rWorkflowVersion;
	private String rCreateTime;
	private String rLastmodTime;
	private String rProgress;
	private Time   rLastModTime;
	// TODO (PDE-577) we need a boolean that would flag the state of a Report
	// (set it to true depending on rLastmodTime)
	// Later we can use this boolean flag for highlighting updated Reports
	private boolean rUpSinceLastTime;
	
	public static final String EMPTY_REPORT = "Empty Report";

	public Report(String sname, String wname, String wversion, String ctime,
			String ltime, boolean updated) {
		rSampleName = sname;
		rWorkflowName = wname;
		rWorkflowVersion = wversion;
		rCreateTime = ctime;
		rLastmodTime = ltime;
		rProgress = null;
		rUpSinceLastTime = updated;
		setTimeStamp(new Time());
		
		if (rWorkflowName!= EMPTY_REPORT){
			String lmTime = rLastmodTime.replaceAll("-", ":").replaceAll(":", "");
			try {
				String parsable = lmTime.substring(0, lmTime.lastIndexOf("."))
						.replace(" ", "T");
				getTimeStamp().parse(parsable);
			} catch (TimeFormatException tfe) {
				getTimeStamp().setToNow();
				Log.e(ReporterActivity.TAG, "An error with Parsing Report Time occured");
			}
		}
	}

	public String getrSampleName() {
		return rSampleName;
	}

	public String getrLastmodTime() {
		return rLastmodTime;
	}

	public String getrWorkflowName() {
		return rWorkflowName;
	}

	public String getrWorkflowVersion() {
		return rWorkflowVersion;
	}

	public String getrCreateTime() {
		return rCreateTime;
	}

	public void setrProgress(String progress) {
		rProgress = progress;
	}

	public String getrProgress() {
		return rProgress;
	}

	public boolean getrUpdated() {
		return rUpSinceLastTime;
	}

	public int progressValue() {
		int progress = 0;
		if (null == rProgress)
			return progress;

		try {
			progress = Integer.parseInt(rProgress);
		} catch (NumberFormatException ne) {
			ne.printStackTrace();
		}
		return progress;
	}
	
	public String toString() {
		return "Report for " + this.getrSampleName();
	}

	public Time getTimeStamp() {
		return rLastModTime;
	}

	private void setTimeStamp(Time rLastModTime) {
		this.rLastModTime = rLastModTime;
	}
	
	
	

}
