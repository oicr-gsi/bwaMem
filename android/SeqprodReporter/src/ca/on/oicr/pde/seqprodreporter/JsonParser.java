package ca.on.oicr.pde.seqprodreporter;

import java.util.ArrayList;
import java.util.List;

import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;
import org.json.JSONTokener;

import android.text.format.Time;
import android.util.Log;
import android.util.TimeFormatException;

/**
 * JsonParser class is for parsing JSON strings
 * @author pruzanov
 *
 */
public class JsonParser {

	private Time lastUpdateTime;
	private Time failedItemUpdateTime;
	private List<Report> parsedJSON;

	public JsonParser(String JsonString, String[] types, Time lastUpdate) {
		this.parsedJSON = new ArrayList<Report>();
		this.lastUpdateTime = null == lastUpdate ? null : new Time();
		this.failedItemUpdateTime = new Time();
		try {
			JSONObject object = (JSONObject) new JSONTokener(JsonString)
					.nextValue();
			for (int t = 0; t < types.length; t++) {
			JSONArray jsonRecords = object.getJSONArray(types[t]);
			for (int i = 0; i < jsonRecords.length(); i++) {
				JSONObject tmp = (JSONObject) jsonRecords.get(i);

				boolean updated = false;

				// 2014-03-21 14:32:23.729
				try {
					Time recordTime = new Time();
					recordTime.set(tmp.getLong("lmtime"));
					if (null != lastUpdate
							&& recordTime.after(lastUpdate))
						updated = true;
					if (null == this.lastUpdateTime || this.lastUpdateTime.before(recordTime))
						this.lastUpdateTime = recordTime;
					
					if (types[t].equals(ReporterActivity.types[1]) 
							&& this.failedItemUpdateTime.before(recordTime)){
						this.failedItemUpdateTime = recordTime;
					}
				} catch (TimeFormatException tfe) {
					Log.e(ReporterActivity.TAG,
							"An error with Parsing Time occured");
				}

				Report newReport = new Report(tmp.getString("sample"),
											  tmp.getString("workflow"),
											  tmp.getString("version"),
											  tmp.getLong("crtime"),
											  tmp.getLong("lmtime"),
											  tmp.getString("wrun_id"),
											  tmp.getString("status"),
											  types[t],
											  updated);
				if (types[t].equals(ReporterActivity.types[2])
						&& tmp.has("progress")) {
					String p = tmp.getString("progress");
					if (null != p && !p.isEmpty())
						newReport.setrProgress(p);
				}

				this.parsedJSON.add(newReport);
			}
			}
			
		} catch (JSONException e) {
			Log.e(ReporterActivity.TAG, "A JSON parsing error occured");
			e.printStackTrace();
		}
	}

	public Time getNewUpdateTime() {
		return lastUpdateTime;
	}

	public List<Report> getParsedJSON() {
		return parsedJSON;
	}
	
	public Time getFailedItemUpdateTime(){
		return failedItemUpdateTime;
	}

}
