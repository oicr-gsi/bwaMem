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

public class JsonLoaderTask extends AsyncTask<String, Void, List<Report>> {

	private WeakReference<ReportListFragment> mParent;
	private String TYPE;
	private Time lastUpdated;
	
	private static final int WEEK   = 7 * 24 * 3600 * 1000;
	private static final int MONTH  = 30 * 24 * 3600 * 1000;
	private static final int YEAR   = 365 * 24 * 3600 * 1000;
	private static final int DECADE = 10 * 365 * 24 * 3600 * 1000;

	public JsonLoaderTask(ReportListFragment parent, String type,
			Time lastUpdate) {
		super();
		mParent = new WeakReference<ReportListFragment>(parent);
		TYPE = type;
		this.lastUpdated = lastUpdate;
	}

	@Override
	protected List<Report> doInBackground(String... params) {
		List<Report> reports = null;

		try {
			reports = getReportsFromDB(params); 
		} catch (NullPointerException npe) {
			Log.e(ReporterActivity.TAG,"There was an error reading database");
		} catch (RuntimeException rte){
			Log.d(ReporterActivity.TAG, "Caught exception");
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

	private List<Report> getReportsFromDB(String... params)
			throws NullPointerException {
		// TODO PDE-650: Need to change this to enable time-constrained queries
		List<Report> results = this.queryReportData(params[0],params[1]);
		return results;
	}

	/**
	 *  Functions for getting data from cursor
	 */
	private ArrayList<Report> queryReportData(String filterWord, String timeRange) {
		
		Cursor result;
		long earliest =  Long.valueOf(System.currentTimeMillis());
		
		if (null != timeRange && !timeRange.isEmpty()) {	    
		    // Not clear at this point what is the best way to handle this other than
		    // by using long if-else block
		    if (timeRange.equals("decade")) {
		      earliest -= DECADE;		    	
		    } else if (timeRange.equals("year")) {
		      earliest -= YEAR;
		    } else if (timeRange.equals("month")) {
		      earliest -= MONTH;
		    } else { // default to 'week'
		      earliest -= WEEK;
		    }	    
		} else {
			earliest -= WEEK;
		}
		
		if (null != filterWord && !filterWord.isEmpty()) {
		 filterWord = "%" + filterWord + "%";
		 result = mParent.get().getActivity().getApplication()
				  .getContentResolver()
				  .query(DataContract.CONTENT_URI, null, 
						 DataContract.WR_TYPE + "= ? AND " + 
				         DataContract.LM_TIME + "> ? " + " AND (" + 
						 DataContract.SAMPLE + " LIKE ? OR " + 
				         DataContract.WORKFLOW + " LIKE ? )", new String[]{TYPE, "" + earliest, filterWord, filterWord}, null);
		} else {
		  result = mParent.get().getActivity().getApplication()
				  .getContentResolver()
				  .query(DataContract.CONTENT_URI, null, DataContract.WR_TYPE + "=?", new String[]{TYPE}, null);
	    }
		ArrayList<Report> rValue = new ArrayList<Report>();

		if (result != null) {
			if (result.moveToFirst()) {
				Time newLatest = null;
				do {
					Report newEntry = getReportDataFromCursor(result);
					Time entryLMTime = newEntry.getTimeStamp();
					if (null != this.lastUpdated 
							&& entryLMTime.after(this.lastUpdated)){
						newEntry.setrUpdated(true); 
					}
					rValue.add(newEntry);
					if (null == newLatest || newLatest.before(entryLMTime)){
						newLatest = entryLMTime;
					}
				} while (result.moveToNext() == true);
				//Initially update the fragment's update times
				if (null == lastUpdated){
					this.mParent.get().setLastUpdateTime(newLatest);
				} 
				//Update the lastUpdated for parent after all items are created
				else if (this.lastUpdated.before(newLatest)){
					// TODO: this code was used here for some reason, need to see if it's abscence does anything
					// if (this.mParent.get().getSectionNumber()-1
					// 8== this.mParent.get().getActivity().getActionBar().getSelectedNavigationIndex()){					
						this.mParent.get().setLastUpdateTime(newLatest);
				}
			
			result.close();
			}
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
		// TODO PDE-650 change String to long
		//String ctime = cursor.getString(cursor
		//		.getColumnIndex(DataContract.CR_TIME));
		//String ltime = cursor.getString(cursor
		//		.getColumnIndex(DataContract.LM_TIME));
		long ctime = cursor.getLong(cursor
				.getColumnIndex(DataContract.CR_TIME));
		long ltime = cursor.getLong(cursor
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


	/**
	 * getReportsFromFile was initially used for loading data from a static json file,
	 * is still here just in case we need some debugging enabled
	 * 
	 * @param params
	 * @return
	 * @throws IOException
	 */
	@SuppressWarnings("unused")
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
		// Insert data into DBsavedInstanceState
		if (null != results) {
			ContentResolver cr = mParent.get().getActivity().getApplication()
					.getContentResolver();
			for (Report rep : results) {
				cr.insert(DataContract.CONTENT_URI, Report.convertToCV(rep));
			}
		}

		return results;
	}

    /**
     * Get list of Report objects from a long JSON string
     */
	@Deprecated
	public List<Report> getRecordsFromJSON(String JsonString, String type) {

		String [] types = {type};
		JsonParser jp = new JsonParser(JsonString, types, this.lastUpdated);
		List<Report> result = jp.getParsedJSON();
		Time newLatest = jp.getNewUpdateTime();

		if (null == lastUpdated){
			this.mParent.get().setLastUpdateTime(newLatest);
		}
		else if (lastUpdated.before(newLatest)){
			
			// Called only when the corresponding fragment's tab is selected
			if (this.mParent.get().getSectionNumber() - 1 == this.mParent.get()
					.getActivity().getActionBar().getSelectedNavigationIndex()) {

				this.mParent.get().setLastUpdateTime(newLatest);
				Log.d(ReporterActivity.TAG, "Updated last update time for "	+ TYPE);
			}
		}

		return result;
	}

}
