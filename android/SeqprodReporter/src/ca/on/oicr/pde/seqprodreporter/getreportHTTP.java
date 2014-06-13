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
	private Context mContext;
	private AndroidHttpClient mClient;
	private static final String SCRIPT = "/getReport.pl?range=";

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

			byte[] byteArray = dbResponse.getBytes();
			String JsonString = new String(byteArray, "UTF-8");
			JsonParser jp = new JsonParser(JsonString, ReporterActivity.types, null);
			List <Report> results = jp.getParsedJSON();
			SharedPreferences sp = mContext.getSharedPreferences(ReporterActivity.PREFERENCE_FILE, Context.MODE_PRIVATE);
			sp.edit().putString("updateTime", ReporterActivity.timeToStringConverter(jp.getNewUpdateTime())).apply();
			//Insert data into db
			//int count = 0; // WE MAY NEED THIS TO DEBUG FURTHER
			if (null != results) {
			    ContentResolver cr = mContext.getApplicationContext().getContentResolver();
				for (Report rep : results) {
					//DEBUG
					//if (rep.getrWorkflowRunType().equals(ReporterActivity.types[2]))
					//	Log.d(ReporterActivity.TAG, "Inserting pending report with progress %" + rep.getrProgress());
					cr.insert(DataContract.CONTENT_URI, Report.convertToCV(rep));
					//count++;
				}
				//Log.d(ReporterActivity.TAG,results.size() + " Reports received from server, " + count + "insert operations");
				return Boolean.TRUE;
			}

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
