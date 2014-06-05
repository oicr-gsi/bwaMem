package ca.on.oicr.pde.seqprodreporter;

import ca.on.oicr.pde.seqprodprovider.DataContract;
import android.content.ContentValues;
import android.text.format.Time;
import android.util.Log;
import android.util.TimeFormatException;

public class Report {
	//TODO need to add status
	private String rSampleName; 
	private String rSequencerRunName;
	private String rStudyName;
	private String rSeqwareAccession;
	private String rWorkflowName;
	private String rWorkflowVersion;
	private String rCreateTime;
	private String rLastmodTime;
	private String rProgress;
	private Time   rLastModTime;
	private boolean rUpSinceLastTime;
	private int _ID;
	
	/*
	 *  This is MySQL - specific, will set this when inserting into database
	 */
	public int getID() {
		return _ID;
	}

	public void setsID(int _ID) {
		this._ID = _ID;
	}

	public static final String EMPTY_REPORT = "Empty Report";
	
	//TODO Need to update constructor once back end of app is modified
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
		
		if (!rWorkflowName.equals(EMPTY_REPORT)){
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
	
	public Report (ContentValues value) {
		rSampleName = value.getAsString(DataContract.SAMPLE);
		rWorkflowName = value.getAsString(DataContract.WORKFLOW);
		rWorkflowVersion = value.getAsString(DataContract.WF_VERSION);
		rCreateTime = value.getAsString(DataContract.CR_TIME);
		rLastmodTime = value.getAsString(DataContract.LM_TIME);
		rProgress = value.getAsString(DataContract.PROGRESS);
		setTimeStamp(new Time());
		
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
	
	public String getrSequencerRunName(){
		return rSequencerRunName;
	}
	
	public String getrStudyName(){
		return rStudyName;
	}
	
	public String getrSeqwareAccession(){
		return rSeqwareAccession;
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
