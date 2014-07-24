var mongoose = require('mongoose');
var Schema = mongoose.Schema;

var reportSchema = new Schema({
	sample_name : String,
	workflow_name : String,
	workflow_version : String,//Number,
	workflow_run_id : String,//Number,
	workflow_run_type : String,
	workflow_run_status : String,
	create_time : String,//Date,
	last_modified_time : String//Date,
});

module.exports = mongoose.model('Report', reportSchema);
