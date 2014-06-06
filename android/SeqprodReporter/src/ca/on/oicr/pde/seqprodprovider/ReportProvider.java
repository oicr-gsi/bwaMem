package ca.on.oicr.pde.seqprodprovider;

import java.util.HashMap;
import java.util.Map;

import android.content.ContentProvider;
import android.content.ContentValues;
import android.content.Context;
import android.database.Cursor;
import android.database.sqlite.SQLiteDatabase;
import android.database.sqlite.SQLiteOpenHelper;
import android.database.sqlite.SQLiteQueryBuilder;
import android.net.Uri;
import android.util.Log;
import ca.on.oicr.pde.seqprodreporter.ReporterActivity;
import ca.on.oicr.pde.seqprodreporter.Report;

public class ReportProvider extends ContentProvider {
	// Helper class for creating/deleting database
	private MainDatabaseHelper mReportHelper;
	private static Map<String, Integer> cachedWorkflowRuns;

	@Override
	public boolean onCreate() {
		mReportHelper = new MainDatabaseHelper(getContext());
		// TODO this may be not such a good place for this container
		cachedWorkflowRuns = new HashMap<String, Integer>();
		return true;
	}

	@Override
	public Cursor query(Uri uri, String[] projection, String selection,
			String[] selectionArgs, String sortOrder) {
		SQLiteQueryBuilder qBuilder = new SQLiteQueryBuilder();

		qBuilder.setTables(DataContract.DATA_TABLE);
		SQLiteDatabase db = this.mReportHelper.getReadableDatabase();
		Cursor c = qBuilder.query(db, projection, selection, selectionArgs,
				null, // groupBy
				null, // having
				sortOrder);
		// TODO need to return cursor with 0 rows if there's no data
		c.setNotificationUri(getContext().getContentResolver(), uri);
		return c;
	}

	@Override
	public synchronized String getType(Uri uri) {
		String contentType = DataContract.CONTENT_ITEM_TYPE;
		if (isTableUri(uri)) {
			contentType = DataContract.CONTENT_DIR_TYPE;
		}
		return contentType;
	}

	@Override
	public synchronized Uri insert(Uri uri, ContentValues value) {
        // This is a 'smart' insert, if it finds that SWID exists already,
		// the record will be updated
		if (null != value && value.containsKey(DataContract.WR_ID)) {
			Report report = new Report(value);
			if (cachedWorkflowRuns.containsKey(report.getrWorkflowRunId())) {
				String[] wrids = new String[1];
				wrids[0] = report.getrWorkflowRunId();
				//Log.d(ReporterActivity.TAG,"Will update sample " + report.getrSampleName());
				update(uri, value, DataContract.WR_ID + "=?", wrids);
			} else {
				SQLiteDatabase db = this.mReportHelper.getWritableDatabase();
				long rowId = db.insert(DataContract.DATA_TABLE, null, value);
				if (rowId > 0) {
					// row inserted, update cached SQW ACCESSION hashmap
					cachedWorkflowRuns.put(report.getrWorkflowRunId(), 1);
					return Uri.withAppendedPath(DataContract.CONTENT_URI,
							String.valueOf(rowId));
				}
			}
		}
		return null;
	}

	@Override
	public synchronized int delete(Uri uri, String selection,
			String[] selectionArgs) {
		SQLiteDatabase db = mReportHelper.getWritableDatabase();
		return db.delete(DataContract.DATA_TABLE, selection, selectionArgs);
	}

	@Override
	public int update(Uri uri, ContentValues values, String selection,
			String[] selectionArgs) {
		// Update Only things that change (lm time, progress, status)
		SQLiteDatabase db = mReportHelper.getWritableDatabase();
		ContentValues selectedValues = new ContentValues();

		selectedValues.put(DataContract.LM_TIME,
				values.getAsString(DataContract.LM_TIME));
		selectedValues.put(DataContract.STATUS,
				values.getAsString(DataContract.STATUS));
		selectedValues.put(DataContract.PROGRESS,
				values.getAsString(DataContract.PROGRESS));

		return db.update(
				DataContract.DATA_TABLE,
				values,
				selection,
				selectionArgs);
	}

	protected static final class MainDatabaseHelper extends SQLiteOpenHelper {

		/*
		 * Instantiates an open helper for the provider's SQLite data repository
		 * Do not do database creation and upgrade here.
		 */
		MainDatabaseHelper(Context context) {
			super(context, DataContract.DATABASE, null, 1);
		}

		/*
		 * Creates the data repository. This is called when the provider
		 * attempts to open the repository and SQLite reports that it doesn't
		 * exist.
		 */
		public void onCreate(SQLiteDatabase db) {

			Log.d(ReporterActivity.TAG, "Database does not exist, Creating " + DataContract.DATABASE);
			db.execSQL(DataContract.REPORT_DB_CREATE);
		}

		@Override
		public void onUpgrade(SQLiteDatabase db, int oldVersion, int newVersion) {
			// Drops than creates main table
			db.execSQL("DROP TABLE IF EXISTS " + DataContract.DATA_TABLE);
			onCreate(db);
		}
	}

	// Does last segment of the Uri match a string of digits?
	/*
	 * private boolean isItemUri(Uri uri) { return
	 * uri.getLastPathSegment().matches("\\d+"); }
	 */

	// Is the last segment of the Uri the name of the data table?
	private boolean isTableUri(Uri uri) {
		return uri.getLastPathSegment().equals(DataContract.DATA_TABLE);
	}
}
