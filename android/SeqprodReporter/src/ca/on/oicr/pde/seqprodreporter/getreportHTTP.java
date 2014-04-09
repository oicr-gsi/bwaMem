package ca.on.oicr.pde.seqprodreporter;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.http.HttpResponse;
import org.apache.http.client.ClientProtocolException;
import org.apache.http.client.ResponseHandler;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.BasicResponseHandler;
import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;
import org.json.JSONTokener;

import android.content.SharedPreferences;
import android.net.http.AndroidHttpClient;
import android.os.AsyncTask;

public class getreportHTTP extends AsyncTask <String, Void, List<Report>> {
	private final String URL;
	private final String Range;
	
	// TODO get these values from shared preferences
	public getreportHTTP(SharedPreferences sp) {
		//TODO append missing portion of the string to URL, if needed
	  	this.URL = sp.getString("pref_hostName", null);
	  	this.Range = sp.getString("prefs_Scope","week");
	}
	
	AndroidHttpClient mClient = AndroidHttpClient.newInstance(ReporterActivity.TAG);
	
	@Override
	protected List<Report> doInBackground(String... params) {
		// TODO Auto-generated method stub
		HttpGet request = new HttpGet(URL);
		JSONResponseHandler responseHandler = new JSONResponseHandler();
		try {
			return mClient.execute(request, responseHandler);
		} catch (ClientProtocolException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return null;
	}
	
	private class JSONResponseHandler implements ResponseHandler<List<Report>> {

		private static final String LONGITUDE_TAG = "lng";
		private static final String LATITUDE_TAG = "lat";
		private static final String MAGNITUDE_TAG = "magnitude";
		private static final String EARTHQUAKE_TAG = "earthquakes";

		@Override
		public List<Report> handleResponse(HttpResponse response)
				throws ClientProtocolException, IOException {
			List<Report> result = new ArrayList<Report>();
			String JSONResponse = new BasicResponseHandler()
					.handleResponse(response);
			try {

				// Get top-level JSON Object - a Map
				JSONObject responseObject = (JSONObject) new JSONTokener(
						JSONResponse).nextValue();

				// Extract value of "earthquakes" key -- a List
				JSONArray earthquakes = responseObject
						.getJSONArray(EARTHQUAKE_TAG);

				// Iterate over earthquakes list
				for (int idx = 0; idx < earthquakes.length(); idx++) {

					// Get single earthquake data - a Map
					JSONObject earthquake = (JSONObject) earthquakes.get(idx);

				}
			} catch (JSONException e) {
				e.printStackTrace();
			}
			return result;
		}
	}


	
	
}
