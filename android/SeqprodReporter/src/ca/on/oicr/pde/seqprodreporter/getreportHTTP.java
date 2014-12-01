package ca.on.oicr.pde.seqprodreporter;

import java.io.IOException;
import java.util.List;

import org.apache.http.HttpResponse;
import org.apache.http.client.ClientProtocolException;
import org.apache.http.client.ResponseHandler;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.BasicResponseHandler;

import android.content.ContentResolver;
import android.content.ContentValues;
import android.content.Context;
import android.content.Intent;
import android.net.http.AndroidHttpClient;
import android.os.AsyncTask;
import android.support.v4.content.LocalBroadcastManager;
import android.text.format.Time;
import android.util.Log;
import ca.on.oicr.pde.seqprodprovider.DataContract;

/**
 * An asynchronous task that is invoked to download the Worklow data into the application's SQLite database
 * from the JSON script. This task is done on a background thread.
 * Once the task is finished it will send a notification to the user indicating that a
 * data update has finished. The form of the notification will vary depending on the state of the App.
 * 
 *@see JsonLoaderTask
 */
public class getreportHTTP extends AsyncTask<Time, Void, Boolean> {
	private final String URL;
	private final String Range;
	private Context mContext;
	private AndroidHttpClient mClient;
	//TODO this needs to be moved to strings.xml
	private static final String SCRIPT = "/getReport.pl?range=";
	
	private Time failedTime;
	private boolean isFailedModified;
	private String updateTime;

	/**
	 * Constructs an instance of getreportHTTP with the application context, the host URL to download
	 * the data from, and the date range to date back to when receiving workflow reports. The workflow
	 * data that is being downloaded will depend on what updateRange value is being passed in. The different
	 * values that could be passed in as the updateRange parameter are: 'week', 'month', 'year' and 'decade'.
	 * For example, if updateRange has the value 'week'; that instance of getreportHTTP will download only
	 * the data of workflows that have been updated within the last week.
	 * @param context the application context of the calling application
	 * @param hostURL the host URL of the JSON script which will be used to download the data from
	 * @param updateRange specifies the range to date back to when getting workflow data
	 */
	public getreportHTTP(Context context, String hostURL, String updateRange) {

		this.URL = hostURL;
		this.Range = updateRange;
		this.mContext = context;
	}
	@Override
	protected Boolean doInBackground(Time... params) {

		if (null == this.URL || null == this.Range)
			return null;
		if (null != params && params.length >= 1)
			failedTime = params[0];
		
		isFailedModified = false;
		Boolean result = Boolean.FALSE;
		String fullURL = this.URL + SCRIPT + this.Range;
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
			Log.d(ReporterActivity.TAG,"Http response received...");
			String dbResponse = new BasicResponseHandler().handleResponse(response);	
			String JsonString = new String(dbResponse.getBytes(), "UTF-8");
			JsonParser jp = new JsonParser(JsonString, ReporterActivity.types, null);
			List <Report> results = jp.getParsedJSON();

			if (null == failedTime || failedTime.before(jp.getFailedItemUpdateTime())){
				failedTime = jp.getFailedItemUpdateTime();
				isFailedModified = true;
			}
			updateTime = ReporterActivity.timeToStringConverter(jp.getNewUpdateTime());
			

			//Insert data into db
			if (null != results) {
				ContentResolver cr = mContext.getApplicationContext().getContentResolver();
				// Removing all entries for 'pending' entries, but...
				// if we don't have any results - don't do it
				String[] pendingWRuns = {ReporterActivity.types[2]};
				int deletedWRs = cr.delete(DataContract.CONTENT_URI,
						                   DataContract.WR_TYPE + "=?",
						                   pendingWRuns);
				if (deletedWRs > 0) {
					Log.d(ReporterActivity.TAG, "Deleted " + deletedWRs + " Old Pending Runs");
				}
			    
				//We use bulk insert instead of incremental insert
			    ContentValues[] dataChunks = new ContentValues[results.size()];
			    for (int r = 0; r < results.size(); r++) {
			      dataChunks[r] = Report.convertToCV(results.get(r));
			    }
			    cr.bulkInsert(DataContract.CONTENT_URI, dataChunks);
				return Boolean.TRUE;
			}

			return Boolean.FALSE;
		}
	}

	@Override
	protected void onPostExecute(Boolean result) {
		Log.d(ReporterActivity.TAG,"Http Request execution finished...");
		Intent intent = new Intent(ReporterActivity.DATACHANGE_INTENT);
		if (result){
			intent.putExtra("updateTime", updateTime);
			if (isFailedModified){
				intent.putExtra("modifiedFailedTime", failedTime.format2445());
			}
		}
		LocalBroadcastManager.getInstance(mContext).sendBroadcast(intent);
	}
	
}
