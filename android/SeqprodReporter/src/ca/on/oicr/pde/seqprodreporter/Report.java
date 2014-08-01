package ca.on.oicr.pde.seqprodreporter;

import ca.on.oicr.pde.seqprodprovider.DataContract;
import android.content.ContentValues;
import android.text.format.Time;
import android.util.Log;
import android.util.TimeFormatException;

public class Report {

	private String rSampleName;
	private String rSequencerRunName;
	private String rStudyName;
	private String rSeqwareAccession;
	private String rWorkflowName;
	private String rWorkflowVersion;
	private String rWorkflowRunId;
	private String rWorkflowRunType; // get this from web in json format
	private String rWorkflowRunStatus; // what we get from db (status field)
	private String rCreateTime;
	private String rLastmodTime;
	private String rProgress;
	private Time rLastModTime;
	private boolean rUpSinceLastTime;
	private int _ID;

	/*
	 * This _ID thing is MySQLite - specific, will set this when inserting into
	 * database
	 */
	public int getID() {
		return _ID;
	}

	public void setsID(int _ID) {
		this._ID = _ID;
	}

	public static final String EMPTY_REPORT = "Empty Report";
	public static final String ZERO_PROGRESS = "0";

	// TODO place to change if adding new fields
	public Report(String sname, String wname, String wversion, String ctime,
			String ltime, String wrun_id, String wrun_status, String wrun_type,
			boolean updated) {
		rSampleName = sname;
		rWorkflowName = wname;
		rWorkflowVersion = wversion;
		rWorkflowRunId = wrun_id;
		rWorkflowRunType = wrun_type;
		rWorkflowRunStatus = wrun_status;
		rCreateTime = ctime;
		rLastmodTime = ltime;
		rProgress = null;
		rUpSinceLastTime = updated;
		setTimeStamp(new Time());

		//TODO need to consolidate code and make this conversion a function?
		if (!rWorkflowName.equals(EMPTY_REPORT)) {
			String lmTime = rLastmodTime.replaceAll("-", ":").replaceAll(":",
					"");
			try {
				String parsable = lmTime.substring(0, lmTime.lastIndexOf("."))
						.replace(" ", "T");
				getTimeStamp().parse(parsable);
			} catch (TimeFormatException tfe) {
				getTimeStamp().setToNow();
				Log.e(ReporterActivity.TAG,
						"An error with Parsing Report Time occured");
			}
		}
	}

	public Report(ContentValues value) {
		rSampleName = value.getAsString(DataContract.SAMPLE);
		rWorkflowName = value.getAsString(DataContract.WORKFLOW);
		rWorkflowVersion = value.getAsString(DataContract.WF_VERSION);
		rWorkflowRunId = value.getAsString(DataContract.WR_ID);
		rWorkflowRunType = value.getAsString(DataContract.WR_TYPE);
		rWorkflowRunStatus = value.getAsString(DataContract.STATUS);
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
			Log.e(ReporterActivity.TAG,
					"An error with Parsing Report Time occured");
		}
	}

	public void setrUpdated(boolean bool){
		rUpSinceLastTime = bool;
	}
	
	public String getrSequencerRunName() {
		return rSequencerRunName;
	}

	public String getrStudyName() {
		return rStudyName;
	}

	public String getrSeqwareAccession() {
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

	public String getrWorkflowRunId() {
		return rWorkflowRunId;
	}

	public String getrWorkflowRunStatus() {
		return rWorkflowRunStatus;
	}

	public String getrWorkflowRunType() {
		return rWorkflowRunType;
	}

	public String getrCreateTime() {
		return rCreateTime;
	}

	public void setrProgress(String progress) {
		if (null == progress)
			rProgress = ZERO_PROGRESS;
		else
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

		try {
			progress = Integer.parseInt(rProgress);
		} catch (NumberFormatException ne) {
			Log.e(ReporterActivity.TAG, "Invalid value for progress detected");
		}
		return progress;
	}

	public String toString() {
		return "Report for " + this.getrSampleName() + " ["
				+ this.getrWorkflowRunStatus() + "]";
	}

	public Time getTimeStamp() {
		return rLastModTime;
	}

	private void setTimeStamp(Time rLastModTime) {
		this.rLastModTime = rLastModTime;
	}

	public static ContentValues convertToCV(Report report) {
		if (null != report) {
			ContentValues values = new ContentValues();
			values.put(DataContract.SAMPLE, report.getrSampleName());
			values.put(DataContract.WORKFLOW, report.getrWorkflowName());
			values.put(DataContract.WF_VERSION, report.getrWorkflowVersion());
			values.put(DataContract.WR_ID, report.getrWorkflowRunId());
			values.put(DataContract.WR_TYPE, report.getrWorkflowRunType());
			values.put(DataContract.STATUS, report.getrWorkflowRunStatus());
			values.put(DataContract.CR_TIME, report.getrCreateTime());
			values.put(DataContract.LM_TIME, report.getrLastmodTime());

			String progress = report.getrProgress();
			if (null != progress) {
				values.put(DataContract.PROGRESS, progress);
			} else {
				values.put(DataContract.PROGRESS, ZERO_PROGRESS);
			}
			return values;
		}
		return null;
	}

}
