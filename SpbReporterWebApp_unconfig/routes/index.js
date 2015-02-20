var express = require('express');
var router = express.Router();
var Report = require('../models/reports');
var async = require('async');
/* GET home page. */
router.get('/', function(req, res) {
	// We will first get the completed, then failed and pending reports from the database each in array form.
	// Then we will pass to the calling function a 2 dimensional array 'reportsArray' which will have the queried results
	// in the first 3 indexes.  Finally we will render the 'index' html file (found in views/index.jade) while assigning
	// variables in the jade file to represent each array that has the queried reports
	async.series([getCompletedReports,getFailedReports,getPendingReports],
		function(err, reportsArray){
			res.render('index',{completed_reports:reportsArray[0], failed_reports:reportsArray[1], pending_reports:reportsArray[2], JSONdownloadTime:downloadTime});
		});

});
module.exports = router;

function getCompletedReports(callback){
	Report.find({workflow_run_type:'completed'}).sort({last_modified_date : 'desc'}).exec(function(err, completed_reports){
		console.log('completed: ' + completed_reports.length);

		if (err){
			callback(err, null);
		}
		else {
			callback(null, completed_reports);
		}
	});
}

function getFailedReports(callback){
	Report.find({workflow_run_type:'failed'}).sort({last_modified_date : 'desc'}).exec(function(err, failed_reports){
		console.log('failed: ' + failed_reports.length);
		if (err){
			callback(err, null);
		}
		else {
			callback(null, failed_reports);
		}
	});
}

function getPendingReports(callback){
	Report.find({workflow_run_type:'pending'}).sort({last_modified_date : 'desc'}).exec(function(err, pending_reports){
		console.log('pending: ' + pending_reports.length);
		if (err){
			callback(err, null);
		}
		else {
			callback(null, pending_reports);
		}
	});
}
