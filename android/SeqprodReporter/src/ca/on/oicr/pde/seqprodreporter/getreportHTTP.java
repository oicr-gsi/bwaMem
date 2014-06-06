package ca.on.oicr.pde.seqprodreporter;

import java.io.IOException;
import java.util.List;

import org.apache.http.HttpResponse;
import org.apache.http.client.ClientProtocolException;
import org.apache.http.client.ResponseHandler;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.BasicResponseHandler;

import android.content.ContentResolver;
import android.content.Context;
import android.content.Intent;
import android.content.SharedPreferences;
import android.net.http.AndroidHttpClient;
import android.os.AsyncTask;
import android.support.v4.content.LocalBroadcastManager;
import ca.on.oicr.pde.seqprodprovider.DataContract;

public class getreportHTTP extends AsyncTask<Void, Void, Boolean> {
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
	protected Boolean doInBackground(Void... params) {

		Boolean result = Boolean.FALSE;
		String fullURL = URL + SCRIPT + this.Range;
		this.mClient = AndroidHttpClient
				.newInstance(ReporterActivity.TAG);
		
		HttpGet request = new HttpGet(fullURL);

		DBResponseHandler responseHandler = new DBResponseHandler();
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

	private class DBResponseHandler implements ResponseHandler<Boolean> {
		
		@Override
		public Boolean handleResponse(HttpResponse response)
				throws ClientProtocolException, IOException {
			String dbResponse = new BasicResponseHandler()
					.handleResponse(response);
			//TODO Here we need to update database, no files anymore
			//ContentResolver cr =  mContext.getApplicationContext().getContentResolver();
			//String FILENAME = ReporterActivity.DATA_FILE
				//	.replace("RANGE", Range);

			byte[] byteArray = dbResponse.getBytes();
			String JsonString = new String(byteArray, "UTF-8");
			JsonParser jp = new JsonParser(JsonString, ReporterActivity.types, null);
			List <Report> results = jp.getParsedJSON();
			//Insert data into db
			if (null != results) {
			    ContentResolver cr = mContext.getApplicationContext().getContentResolver();
				for (Report rep : results) {
					cr.insert(DataContract.CONTENT_URI, Report.convertToCV(rep));
				}
				return Boolean.TRUE;
			}
				/*if (null != records && records.length > 0) {
					FileOutputStream fos = mContext.openFileOutput(FILENAME,
							Context.MODE_PRIVATE);
					fos.write(records);
					fos.close();
					Log.d(ReporterActivity.TAG, "Saved data to File " + FILENAME);
					return Boolean.TRUE;*/
				//}
				//return Boolean.FALSE;


			return Boolean.FALSE;
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
