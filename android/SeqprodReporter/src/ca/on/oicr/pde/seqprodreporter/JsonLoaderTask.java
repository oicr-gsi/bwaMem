package ca.on.oicr.pde.seqprodreporter;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.lang.ref.WeakReference;
import java.util.ArrayList;
import java.util.List;

import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;
import org.json.JSONTokener;

import android.content.res.Resources;
import android.os.AsyncTask;
import android.text.format.Time;
import android.util.Log;
import android.util.TimeFormatException;

public class JsonLoaderTask extends AsyncTask<Void, Void, List<Report>> {

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
	protected List<Report> doInBackground(Void... params) {
		List<Report> report = null;
		try {
			report = getReportsFromFile();
		} catch (IOException ie) {
			ie.printStackTrace();
		}
		return report;

	}

	@Override
	protected void onPostExecute(List<Report> result) {

		if (null != result && null != mParent.get()) {
			Log.d(ReporterActivity.TAG, "Loaded " + result.size()
					+ " records for " + TYPE);
			mParent.get().addLocalReports(result);
		}
	}

	private List<Report> getReportsFromFile(Void... params) throws IOException {

		Resources r = mParent.get().getResources();
		// TODO TEMPORARY - read only from week' file, implement other time ranges
		String jsonLine = "";
		String jString = "";
		String FNAME = ReporterActivity.DATA_FILE.replace("RANGE", "week");

		boolean localFileOk = false;
		BufferedReader buf = null;

		if (mParent.get().getActivity().getApplicationContext()
				.getFileStreamPath(FNAME).exists()) {
			// First, try to read from locally stored file (previous update)
			try {
				FileInputStream fis = mParent.get().getActivity()
						.getApplicationContext().openFileInput(FNAME);
				buf = new BufferedReader(new InputStreamReader(fis));
				do {
					jString = jString.concat(jsonLine);
					jsonLine = buf.readLine();
				} while (null != jsonLine && !jsonLine.isEmpty());
				buf.close();
				localFileOk = !jString.isEmpty();
			} catch (FileNotFoundException fnfe) {
				Log.d(ReporterActivity.TAG, "Could not find cached data in "
						+ FNAME);
			} finally {
				if (null != buf)
					buf.close();
			}
		}

		// If we don't have update data, read from static file (may be really old)
		if (!localFileOk) {
			Log.d(ReporterActivity.TAG,
					"Will Load the data from static test.json");
			InputStream fis = r.openRawResource(R.raw.test);
			BufferedReader br = new BufferedReader(new InputStreamReader(fis));

			while (null != (jsonLine = br.readLine())) {
				jString = jString.concat(jsonLine);
			}
			br.close();
		} else {
			Log.d(ReporterActivity.TAG, "Loaded cached data");
		}

		return getRecordsFromJSON(jString, this.TYPE);
	}

	public List<Report> getRecordsFromJSON(String JsonString, String type) {
		List<Report> result = new ArrayList<Report>();
		try {
			JSONObject object = (JSONObject) new JSONTokener(JsonString)
					.nextValue();
			JSONArray Reports = object.getJSONArray(type);
			Time newLatest = new Time();
			for (int i = 0; i < Reports.length(); i++) {
				JSONObject tmp = (JSONObject) Reports.get(i);

				String lmTime = tmp.getString("lmtime").replaceAll("-", ":")
						.replaceAll(":", "");
				boolean updated = false;

				// 2014-03-21 14:32:23.729
				try {
					String parsable = lmTime.substring(0,
							lmTime.lastIndexOf(".")).replace(" ", "T");
					Time recordTime = new Time();
					recordTime.parse(parsable);
					if (null != this.lastUpdated
							&& recordTime.after(this.lastUpdated))
						updated = true;
					if (null == newLatest || newLatest.before(recordTime))
						newLatest = recordTime;

				} catch (TimeFormatException tfe) {
					Log.e(ReporterActivity.TAG,
							"An error with Parsing Time occured");
				}

				Report newReport = new Report(tmp.getString("sample"),
						tmp.getString("workflow"), tmp.getString("version"),
						tmp.getString("crtime"), tmp.getString("lmtime"),
						updated);
				if (type.equals(ReporterActivity.types[2])
						&& tmp.has("progress")) {
					String p = tmp.getString("progress");
					if (null != p && !p.isEmpty())
						newReport.setrProgress(p);
				}

				result.add(newReport);
			}
			if (null == lastUpdated || lastUpdated.before(newLatest))

				// Call only when the corresponding fragment is visible
				if (this.mParent.get().isVisible())
					this.mParent.get().setLastUpdateTime(newLatest);

		} catch (JSONException e) {
			e.printStackTrace();
		}
		return result;
	}

}
