var express = require('express');
var router = express.Router();
var Report = require('../models/reports');
var async = require('async');
var workflow_run_types = ["completed","failed", "pending"];

router.get('/', function(req, res) {
	var workflow_stats_data = [];
	Report.distinct('workflow_name',{},function(err, results){
		console.log(results);
		async.each(results,
			function(workflow, callback){
				async.parallel([
					function(callback){
						Report.count({workflow_run_type:'completed', workflow_name:workflow}, function(err, count){
							callback(err, count);
						});
					},
					function(callback){
						Report.count({workflow_run_type:'failed', workflow_name:workflow}, function(err, count){
							callback(err, count);
						});
					},
					function(callback){
						Report.count({workflow_run_type:'pending', workflow_name:workflow}, function(err, count){
							callback(err, count);
						});
					}],
				function(err, workflow_stats){
					if (err){
						console.error(err);
					}
					else {
						workflow_stats_data.push({workflow:workflow, stats:workflow_stats});
					}
					callback(null);
				});
			},
			function(err){
				async.parallel([
					function(callback){
						Report.count({workflow_run_type:'completed'}, function(err, count){
							callback(err,count);
						});
					},
					function(callback){
						Report.count({workflow_run_type:'failed'}, function(err, count){
							callback(err,count);
						});
					},
					function(callback){
						Report.count({workflow_run_type:'pending'}, function(err, count){
							callback(err,count);
						});
					}],
					function(err,all_workflow_stats){
						res.render('stats',{stats:workflow_stats_data, all:all_workflow_stats});
				});
			});
	});
});
module.exports = router;
