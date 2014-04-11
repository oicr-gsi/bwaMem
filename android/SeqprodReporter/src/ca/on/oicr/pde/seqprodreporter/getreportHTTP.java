package ca.on.oicr.pde.seqprodreporter;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;

import org.apache.http.HttpResponse;
import org.apache.http.client.ClientProtocolException;
import org.apache.http.client.ResponseHandler;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.BasicResponseHandler;

import android.content.Context;
import android.content.Intent;
import android.content.SharedPreferences;
import android.net.http.AndroidHttpClient;
import android.os.AsyncTask;
import android.support.v4.content.LocalBroadcastManager;
import android.util.Log;

public class getreportHTTP extends AsyncTask <String, Void, Boolean> {
	private final String URL;
	private final String Range;
	private static final String SCRIPT = "/getReport.pl?range=";
	private Context mContext;
	AndroidHttpClient mClient = AndroidHttpClient.newInstance(ReporterActivity.TAG);
	
	public getreportHTTP(Context context, SharedPreferences sp) {

		this.URL = sp.getString("pref_hostName", null);
	  	this.Range = sp.getString("prefs_Scope","week");
	  	this.mContext = context;
	}
	
	
	
	@Override
	protected Boolean doInBackground(String... params) {

		Boolean result = Boolean.valueOf(false);
		String fullURL = URL + SCRIPT + this.Range;
		HttpGet request = new HttpGet(fullURL);
		JSONResponseHandler responseHandler = new JSONResponseHandler();
		try {
			return mClient.execute(request, responseHandler);
		} catch (ClientProtocolException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (NullPointerException e) {
			e.printStackTrace();
		}
		return result;
	}
	
	private class JSONResponseHandler implements ResponseHandler<Boolean> {

		@Override
		public Boolean handleResponse(HttpResponse response)
				throws ClientProtocolException, IOException {
			Boolean result = Boolean.valueOf(false);
			String JSONResponse = new BasicResponseHandler()
					.handleResponse(response);
			String FILENAME = ReporterActivity.DATA_FILE.replace("RANGE", Range);
			try {
				FileOutputStream fos =  mContext.openFileOutput(FILENAME, Context.MODE_PRIVATE);
				fos.write(JSONResponse.getBytes());
				fos.close();			
				Log.d(ReporterActivity.TAG, "Saved data to File " + FILENAME);
				return Boolean.TRUE; 
				
			} catch (FileNotFoundException e) {
				Log.e(ReporterActivity.TAG, "Could not save to File " + FILENAME);
			}
			return result;
		}
	}
	
	@Override
	protected void onPostExecute(Boolean result) {
        // TODO if result is true, send a Broadcast to notify ReporterActivity
		LocalBroadcastManager.getInstance(mContext)
		 .sendBroadcast(new Intent(ReporterActivity.DATACHANGE_INTENT));
	}
	
}
