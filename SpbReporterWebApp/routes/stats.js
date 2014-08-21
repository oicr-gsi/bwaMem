var express = require('express');
var router = express.Router();
var Report = require('../models/reports');
var async = require('async');
var workflow_run_types = ["completed","failed", "pending"];
// This is the route that 'gets' our stats page
router.get('/', function(req, res) {
	// JSON object that will be used to store the stats for all workflows
	var workflow_stats_data = [];

	// We will first get the distinct workflows that we have in our DB and put each of the names in the 'results' array
	Report.distinct('workflow_name',{},function(err, results){
		console.log(results);
		// For each of the distinct workflows that we have, we will get the number of 'completed', 'failed' and 'pending' that
		// correspond to this workflow in parallel.
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
					//'workflow_stats' is an array that has the number of 'completed','failed' and 'pending' workflows for that specfic workflow
				function(err, workflow_stats){
					if (err){
						console.error(err);
					}
					else {
						// Here we will add onto the JSON with the workflow key storing the workflow name, and the stats key storing the workflow_stats array
						// that contains the stats for that specfic workflow
						workflow_stats_data.push({workflow:workflow, stats:workflow_stats});
					}
					callback(null);
				});
			},
			// Once we are done getting the stats for all the different workflows and placing it in the JSON object, we will
			// get the stats for all workflows using the same method and then finally rendering the stats page (views/stats.jade) while passing in the
			// stats for all workflows as well as all stats to the view
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
