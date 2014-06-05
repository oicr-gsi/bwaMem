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

public class getreportHTTP extends AsyncTask<String, Void, Boolean> {
	private final String URL;
	private final String Range;
	private static final String SCRIPT = "/getReport.pl?range=";
	private Context mContext;
	private AndroidHttpClient mClient;

	public getreportHTTP(Context context, SharedPreferences sp) {

		this.URL = sp.getString("pref_hostName", null);
		this.Range = sp.getString("prefs_Scope", "week");
		this.mContext = context;
	}

	@Override
	protected Boolean doInBackground(String... params) {

		Boolean result = Boolean.valueOf(false);
		String fullURL = URL + SCRIPT + this.Range;
		this.mClient = AndroidHttpClient
				.newInstance(ReporterActivity.TAG);
		
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
		} finally {
		  mClient.close();
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
			//TODO Here we need to update database, no files anymore
			//ContentResolver cr =  mContext.getApplicationContext().getContentResolver();
			String FILENAME = ReporterActivity.DATA_FILE
					.replace("RANGE", Range);
			try {
				byte[] records = JSONResponse.getBytes();
				if (null != records && records.length > 0) {
					FileOutputStream fos = mContext.openFileOutput(FILENAME,
							Context.MODE_PRIVATE);
					fos.write(records);
					fos.close();
					Log.d(ReporterActivity.TAG, "Saved data to File " + FILENAME);
					return Boolean.TRUE;
				}
				return Boolean.FALSE;

			} catch (FileNotFoundException e) {
				Log.e(ReporterActivity.TAG, "Could not save to File "
						+ FILENAME);
			}
			return result;
		}
	}

	@Override
	protected void onPostExecute(Boolean result) {
		// If result is true, send a Broadcast to notify ReporterActivity
		if (result)
		   LocalBroadcastManager.getInstance(mContext).sendBroadcast(
				new Intent(ReporterActivity.DATACHANGE_INTENT));
	}

}
