package ca.on.oicr.pde.seqprodprovider;

import android.content.ContentResolver;
import android.net.Uri;

public final class DataContract {

	public static final String SAMPLE = "sample_name";
	public static final String SEQRUN_NAME = "seqrun_name";
	public static final String STUDY_NAME = "study_name";
	public static final String WR_ID = "workflow_run_id";
	public static final String WR_TYPE = "workflow_run_type";
	public static final String WORKFLOW = "workflow_name";
	public static final String WF_VERSION = "workflow_version";
	public static final String STATUS  = "status";
	public static final String CR_TIME = "create_time";
	public static final String LM_TIME = "lm_time";
	public static final String PROGRESS = "progress";
	// DB parameters
	public static final String DATA_TABLE = "report_table";
	public static final String DATABASE = "reports";
	
	private static final Uri BASE_URI = Uri
			.parse("content://ca.on.oicr.pde.seqprodprovider/");

	// The URI for this table.
	public static final Uri CONTENT_URI = Uri.withAppendedPath(BASE_URI,
			DATA_TABLE);

	// Mime type for a directory of data items
	public static final String CONTENT_DIR_TYPE = ContentResolver.CURSOR_DIR_BASE_TYPE
			+ "/ReportProvider.data.text";

	// Mime type for a single data item
	public static final String CONTENT_ITEM_TYPE = ContentResolver.CURSOR_ITEM_BASE_TYPE
			+ "/ReportProvider.data.text";

	// All columns of this table
	public static final String[] ALL_COLUMNS = { SAMPLE, SEQRUN_NAME, STUDY_NAME, WR_ID, WR_TYPE, WORKFLOW, WF_VERSION, STATUS, CR_TIME, LM_TIME, PROGRESS };


	// Database initialization Query
	public static final String REPORT_DB_CREATE = "CREATE TABLE " +
			DATA_TABLE +
		    "(" +                           
		    " _ID INTEGER PRIMARY KEY AUTOINCREMENT, " + // this is needed to make it all work for ListView
		    SAMPLE + " TEXT, " +
		    SEQRUN_NAME + " TEXT, " +
		    STUDY_NAME + " TEXT," +
		    WR_ID  + " TEXT," +
		    WR_TYPE  + " TEXT, " + // such as pending, completed etc
		    WORKFLOW + " TEXT, " +
		    WF_VERSION + " TEXT, " +
		    STATUS + " TEXT, " +  // actual status (maybe running, scheduled etc.)
		    CR_TIME + " TEXT, " +
		    LM_TIME + " TEXT, " +
 		    PROGRESS + " TEXT, " + 
 		   "UNIQUE (" + DataContract.WR_ID + ") ON CONFLICT REPLACE)";
    }