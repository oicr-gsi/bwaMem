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
import android.net.http.AndroidHttpClient;
import android.os.AsyncTask;
import android.support.v4.content.LocalBroadcastManager;
import android.text.format.Time;
import ca.on.oicr.pde.seqprodprovider.DataContract;

public class getreportHTTP extends AsyncTask<Time, Void, Boolean> {
	private final String URL;
	private final String Range;
	private Context mContext;
	private AndroidHttpClient mClient;
	private static final String SCRIPT = "/getReport.pl?range=";
	
	private Time failedTime;
	private boolean isFailedModified;
	private String updateTime;

	public getreportHTTP(Context context, String hostURL, String updateRange) {

		this.URL = hostURL;
		this.Range = updateRange;
		this.mContext = context;
	}
	@Override
	protected Boolean doInBackground(Time... params) {

		if (null == this.URL || null == this.Range)
			return null;
		if (null != params)
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
			String dbResponse = new BasicResponseHandler()
					.handleResponse(response);
			byte[] byteArray = dbResponse.getBytes();
			String JsonString = new String(byteArray, "UTF-8");
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
				for (Report rep : results) {
					cr.insert(DataContract.CONTENT_URI, Report.convertToCV(rep));
				}
				return Boolean.TRUE;
			}

			return Boolean.FALSE;
		}
	}

	@Override
	protected void onPostExecute(Boolean result) {
		Intent intent = new Intent(ReporterActivity.DATACHANGE_INTENT);
		if (result){
			intent.putExtra("updateTime", updateTime);
			if (isFailedModified){
				intent.putExtra("modifiedFailedTime", failedTime.format2445());
			}
		}
		 LocalBroadcastManager.getInstance(mContext).sendBroadcast(
					intent);
	}
	
}
