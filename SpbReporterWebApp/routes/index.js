var express = require('express');
var router = express.Router();
var Report = require('../models/reports');
var async = require('async');

/* GET home page. */
router.get('/', function(req, res) {
	async.series([getCompletedReports,getFailedReports,getPendingReports ],
		function(err, reportsArray){
			res.render('index',{completed_reports:reportsArray[0], failed_reports:reportsArray[1], pending_reports:reportsArray[2]});
		});

});
module.exports = router;

function getCompletedReports(callback){
	Report.find({workflow_run_type:'completed'}, function(err, completed_reports){
		if (err){
			callback(err, null);
		}
		else {
			callback(null, completed_reports);
		}
	});
}

function getFailedReports(callback){
	Report.find({workflow_run_type:'failed'}, function(err, failed_reports){
		if (err){
			callback(err, null);
		}
		else {
			callback(null, failed_reports);
		}
	});
}

function getPendingReports(callback){
	Report.find({workflow_run_type:'pending'}, function(err, pending_reports){
		if (err){
			callback(err, null);
		}
		else {
			callback(null, pending_reports);
		}
	});
}
