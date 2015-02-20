var mongoose = require('mongoose');
var Schema = mongoose.Schema;

// The schema for our mongoose model that will be exported and used in our App
var reportSchema = new Schema({
	//** Date values are stored in GMT time zone by default**
	sample_name : String,
	workflow_name : String,
	workflow_version : String,
	workflow_run_id : String,
	workflow_run_type : String,
	workflow_run_status : String,
	create_time : String,
	last_modified_time : String,
	progress : String,
	create_date: Date,
	last_modified_date: Date,
	recently_modified: Boolean
});
/**
* This is the model we use to save 'reports' into our MongoDB database, we use this explicitly in our download.js script
* Instances of this model represent a row in our Mongo database.
*/
module.exports = mongoose.model('Report', reportSchema);
