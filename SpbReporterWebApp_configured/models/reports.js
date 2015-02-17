var mongoose = require('mongoose');
var Schema = mongoose.Schema;

var reportSchema = new Schema({
	//** Date values are stored in GMT time zone by default**
	sample_name : String,
	workflow_name : String,
	workflow_version : String,//Number,
	workflow_run_id : String,//Number,
	workflow_run_type : String,
	workflow_run_status : String,
	create_time : String,//Date,
	last_modified_time : String,//Date,
	progress : String,
	create_date: Date,
	last_modified_date: Date,
	recently_modified: Boolean
});
module.exports = mongoose.model('Report', reportSchema);
