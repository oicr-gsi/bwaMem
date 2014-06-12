package ca.on.oicr.pde.seqprodreporter;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.lang.ref.WeakReference;
import java.util.ArrayList;
import java.util.List;

import android.content.ContentResolver;
import android.content.res.Resources;
import android.database.Cursor;
import android.os.AsyncTask;
import android.text.format.Time;
import android.util.Log;
import ca.on.oicr.pde.seqprodprovider.DataContract;

public class JsonLoaderTask extends AsyncTask<Boolean, Void, List<Report>> {

	private WeakReference<ReportListFragment> mParent;
	private String TYPE;
	private Time lastUpdated;

	public JsonLoaderTask(ReportListFragment parent, String type,
			Time lastUpdate) {
		super();
		mParent = new WeakReference<ReportListFragment>(parent);
		TYPE = type;
		this.lastUpdated = lastUpdate;
	}

	@Override
	protected List<Report> doInBackground(Boolean... params) {
		List<Report> reports = null;
		try {
			reports = params[0].booleanValue() == Boolean.TRUE ? getReportsFromFile()
					: getReportsFromDB();
		} catch (NullPointerException npe) { // TODO fix later
			npe.printStackTrace();
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
		return reports;

	}

	@Override
	protected void onPostExecute(List<Report> result) {

		if (null != result && null != mParent.get()) {
			Log.d(ReporterActivity.TAG, "Loaded " + result.size()
					+ " records for " + TYPE);
			mParent.get().addLocalReports(result);
		}
	}

	private List<Report> getReportsFromDB(Void... params)
			throws NullPointerException {
		List<Report> results = this.queryReportData();
		return results;
	}

	/*
	 *  Functions for getting data from cursor
	 */
	private ArrayList<Report> queryReportData() {
		Cursor result = mParent.get().getActivity().getApplication()
				.getContentResolver()
				.query(DataContract.CONTENT_URI, null, DataContract.WR_TYPE + "=?", new String[]{TYPE}, null);
		ArrayList<Report> rValue = new ArrayList<Report>();

		if (result != null) {
			if (result.moveToFirst()) {
				Time newLatest = new Time();
				do {
					Report newEntry = getReportDataFromCursor(result);
					rValue.add(newEntry);
					if (null != this.lastUpdated 
							&& newEntry.getTimeStamp().after(this.lastUpdated)){
						newEntry.setrUpSinceLastTime(true); 
						Log.d(ReporterActivity.TAG, "A report was updated");
					}
					if (null == newLatest || newLatest.before(newEntry.getTimeStamp())){
						newLatest = newEntry.getTimeStamp();
					}
				} while (result.moveToNext() == true);
				
				if (null == lastUpdated || lastUpdated.before(newLatest)){
					
					if (this.mParent.get().getSectionNumber()-1
							== this.mParent.get().getActivity().getActionBar().getSelectedNavigationIndex()){
						
						this.mParent.get().setLastUpdateTime(newLatest);
						Log.d(ReporterActivity.TAG, "Updated last update time for " + TYPE);
					}
				}
			}
			result.close();
		}
		return rValue;
	}

	/*
	 * Create a Report and return it, set progress only for pending runs
	 */
	private Report getReportDataFromCursor(Cursor cursor) {
		String sname = cursor.getString(cursor
				.getColumnIndex(DataContract.SAMPLE));
		String wname = cursor.getString(cursor
				.getColumnIndex(DataContract.WORKFLOW));
		String wversion = cursor.getString(cursor
				.getColumnIndex(DataContract.WF_VERSION));
		String ctime = cursor.getString(cursor
				.getColumnIndex(DataContract.CR_TIME));
		String ltime = cursor.getString(cursor
				.getColumnIndex(DataContract.LM_TIME));
		String wrunid = cursor.getString(cursor
				.getColumnIndex(DataContract.WR_ID));
		String wrunstatus = cursor.getString(cursor
				.getColumnIndex(DataContract.STATUS));
		String wruntype   = cursor.getString(cursor
				.getColumnIndex(DataContract.WR_TYPE));
		String wrprogress = cursor.getString(cursor
				.getColumnIndex(DataContract.PROGRESS));
		// Construct Report from obtained values
		Report newEntry = new Report(sname, wname, wversion, ctime, ltime, wrunid, wrunstatus, wruntype, false);
		if (wruntype.equals(ReporterActivity.types[2]))
			newEntry.setrProgress(wrprogress);
		return newEntry;
	}


	@Deprecated
	private List<Report> getReportsFromFile(Void... params) throws IOException {

		Resources res = mParent.get().getResources();

		String jsonLine = "";
		String jString = "";

		Log.d(ReporterActivity.TAG, "Will Load the data from static test.json");
		InputStream fis = res.openRawResource(R.raw.test);
		BufferedReader br = new BufferedReader(new InputStreamReader(fis));

		while (null != (jsonLine = br.readLine())) {
			jString = jString.concat(jsonLine);
		}
		br.close();
		List<Report> results = getRecordsFromJSON(jString, this.TYPE);
		// Insert data into DB
		if (null != results) {
			ContentResolver cr = mParent.get().getActivity().getApplication()
					.getContentResolver();
			for (Report rep : results) {
				cr.insert(DataContract.CONTENT_URI, Report.convertToCV(rep));
			}
		}

		return results;
	}

    /*
     * Get list of Report objects from a long JSON string
     */
	public List<Report> getRecordsFromJSON(String JsonString, String type) {

		String [] types = {type};
		JsonParser jp = new JsonParser(JsonString, types, this.lastUpdated);
		List<Report> result = jp.getParsedJSON();
		Time newLatest = jp.getNewUpdateTime();

		if (null == lastUpdated || lastUpdated.before(newLatest))
			
			// Called only when the corresponding fragment's tab is selected
			if (this.mParent.get().getSectionNumber() - 1 == this.mParent.get()
					.getActivity().getActionBar().getSelectedNavigationIndex()) {

				this.mParent.get().setLastUpdateTime(newLatest);
				Log.d(ReporterActivity.TAG, "Updated last update time for "	+ TYPE);
			}

		return result;
	}

}
