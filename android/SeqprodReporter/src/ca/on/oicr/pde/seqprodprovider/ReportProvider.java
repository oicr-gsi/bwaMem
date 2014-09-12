package ca.on.oicr.pde.seqprodprovider;

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

public class ReportProvider extends ContentProvider {

	// Helper class for creating/deleting database
	private MainDatabaseHelper mReportHelper;

	@Override
	public boolean onCreate() {
		mReportHelper = new MainDatabaseHelper(getContext());
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
		// TODO need to return cursor with 0 rows if there's no data (although
		// it seems functional regardless)
		c.setNotificationUri(getContext().getContentResolver(), uri);
		return c;
	}

	@Override
	public synchronized int bulkInsert(Uri uri, ContentValues[] values) {
		SQLiteDatabase db = this.mReportHelper.getWritableDatabase();
		db.beginTransactionNonExclusive();
         
        for(int x = 0; x < values.length; x++){
            db.insert(DataContract.DATA_TABLE, null, values[x]);            
        }
 
        db.setTransactionSuccessful();
        db.endTransaction();
        db.close();
		
		return values.length;
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
		if (null != value && value.containsKey(DataContract.WR_ID)) {
			SQLiteDatabase db = this.mReportHelper.getWritableDatabase();
			long rowId = db.insert(DataContract.DATA_TABLE, null, value);
			if (rowId > 0) {
				return Uri.withAppendedPath(DataContract.CONTENT_URI,
						String.valueOf(rowId));
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
	public synchronized int update(Uri uri, ContentValues values, String selection,
			String[] selectionArgs) {
		// Update Only things that change (lm time, progress, status)
		// This method is implemented but NOT IN USE at this point
		SQLiteDatabase db = mReportHelper.getWritableDatabase();
		ContentValues selectedValues = new ContentValues();
		selectedValues.put(DataContract.LM_TIME,
				values.getAsLong(DataContract.LM_TIME));
		selectedValues.put(DataContract.STATUS,
				values.getAsString(DataContract.STATUS));
		selectedValues.put(DataContract.PROGRESS,
				values.getAsString(DataContract.PROGRESS));

		return db.update(DataContract.DATA_TABLE, selectedValues, selection,
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

			Log.d(ReporterActivity.TAG, "Database does not exist, Creating "
					+ DataContract.DATABASE);
			db.execSQL(DataContract.REPORT_DB_CREATE);
		}

		@Override
		public void onUpgrade(SQLiteDatabase db, int oldVersion, int newVersion) {
			// Drops than creates main table
			db.execSQL("DROP TABLE IF EXISTS " + DataContract.DATA_TABLE);
			onCreate(db);
		}
	}

	// Is the last segment of the Uri the name of the data table?
	private boolean isTableUri(Uri uri) {
		return uri.getLastPathSegment().equals(DataContract.DATA_TABLE);
	}
}
