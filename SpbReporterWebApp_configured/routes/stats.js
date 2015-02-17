/* This file retrieves the data acquired from Mongodb. As Mongodb functions need
 * time to proceed and Javascript is synchronous, we need to use async module in JS.
 * The basic idea is wrapping tasks in function objects and executing them in 
 * async.parallel(). When adding new tasks, just push it to big_task[] and use 
 * callbacks to pass data.
 *
 * The structure to represent data of the table is a nested object. Three objects
 * will be generated and sent to the presentation side as table for completed, failed
 * and pending workflows. Three array of short project code will also be rendered in
 * client side to maintain the order of data
 * Example:
 * table
 * {'A' : { 'a' : 1,
 *         'b' : 0 },
 *  'B' : { 'a' : 3,
 *          'b' : 2 }}
 * where 'A', 'B' are workflow names and 'a', 'b' are short project codes
 * short project code array
 * ['a','b']
 */

var express = require('express');
var router = express.Router();
var Report = require('../models/reports');
var async = require('async');

//function that consumes an array of sample_name strings and manipulate the strings.
//remove duplicates and only keep short project code (words before the first underscore)
var cut_strings = function (strl) {
    var newl = [];
    var isIn = function (arr, item) {
        if (arr.length === 0) {
          return false;
        }
        for (var i = 0; i < arr.length; ++i) {
            if (arr[i] === item) {
                return true;
            }
        }
        return false;
    };
    for (var i = 0; i < strl.length; ++i) {
        var cut = strl[i].substring(0, strl[i].indexOf('_'));
        if (!isIn(newl, cut)) {
            newl.push(cut);
        }
    }
    return newl;
};

// This is the route that 'gets' our stats page
router.get('/', function(req, res) {
    var big_task = []; //variable that stores the highest level tasks.
    var workflow_run_types = ["completed", "failed", "pending"];

    //make three function objects which proceed the workflow_name-short_project _code table, corresponding to three workflow_run_types
    //the functions will be stored in big_Task as an array and executed in async.parallel later
    workflow_run_types.forEach(
        function (type) {
            var task_buffer = function(callback2) {
                //For the usage of Report.distinct,see Mongoose API
                Report.distinct('workflow_name', {workflow_run_type : type}, function(err, workflows) {
                    Report.distinct('sample_name', {workflow_run_type : type}, function(err, samples) {
                        var tbuffer = {}; //temporary variable that stores data of the table
                        var sbuffer = []; //temporary variable that stores the short project code

                        //check if there are items
                        if (samples.length > 0 && workflows.length > 0) {
                            console.log(samples);
                            var tasks = []; //store the functions that retrieves data for each combination of workflows and short project codes
                            sbuffer = cut_strings(samples);
                            console.log(sbuffer);

                            //generates the function objects that retrieves data for each combination of workflows and short project codes
                            //by looping through the arrays of them
                            workflows.forEach(function(wf) { // wf is the string of a workflow name
                                tbuffer[wf] = {}; //initialise the property that stores data of wf
                                sbuffer.forEach(
                                    function(sp) { //sp is the string of a short project code
                                        tbuffer[wf][sp] = 0; //initialise the property that stores the count of wf and sp

                                        //RegExp is used to find all workflows whose sample_names begin with sp
                                        var re = new RegExp('^' + sp);

                                        //query is the function object which does the counting job and returns the number
                                        //and corresponding workflow name and short project code as an array of three elements
                                        //to the callback of async.parallel below. wf and sp are to identified the location
                                        //of the number in table
                                        var query = function(callback) {
                                            Report.count({
                                                    workflow_name: wf,
                                                    sample_name: { $regex: re },
                                                    workflow_run_type: type
                                                },
                                                function(err, count) {
                                                    callback(err, [wf, sp, count]);
                                            });
                                        };

                                        //push the function to task array
                                        tasks.push(query);
                                    });
                            });

                            //execute the functions that do the counting job in parallel, the results will be returned to the callback as an array
                            async.parallel(tasks, function(err, count) {

                                //loop that assign the numbers to table buffer
                                for (var i = 0; i < count.length; ++i) {
                                    
				    //temporary variables to simplify nested arrays selector
                                    //the structure of each element in count is
                                    // [workflow_name, short_project_code, actual_count_result]
                                    var tmp = count[i]; //an array of three element
                                    var w = tmp[0]; //workflow name
                                    var s = tmp[1]; //short project code
                                    var c = tmp[2]; //count result
                                    tbuffer[w][s] = c;
                                }

                                //this callback returns table-short_project_codes pair to the final callback of the most outer async.parallel
                                callback2(null, [tbuffer, sbuffer]);
                            });
                        } 
			else { //if there is nothing to retrieve
                            callback2(null, [{}, []]);
                        }
                    });
                });
            };

            //task_buffer is a function object that generates the table of one workflow_run_type
            big_task.push(task_buffer);
    })

    //sum_task is a function object that generates the table of workflow_name and workflow_run_type counts.
    var sum_task = function(callback2) {
        Report.distinct('workflow_name', {}, function(err, results) {
            console.log(results);
            var tasks = [];
            // For each of the distinct workflows that we have, generate the function object that gets the number of 'completed', 'failed' and 'pending' that
            // correspond to this workflow in parallel.
            results.forEach( function(workflow) {
                    var task_buffer = function (callback1) {
                        async.parallel([
                            function(callback) {
                                Report.count({workflow_run_type: 'completed', workflow_name: workflow}, 
				    function(err, count) { callback(err, count); });
                            },
                            function(callback) {
                                Report.count({workflow_run_type: 'failed', workflow_name: workflow}, 
				    function(err, count) { callback(err, count); });
                            },
                            function(callback) {
                                Report.count({workflow_run_type: 'pending', workflow_name: workflow}, 
				    function(err, count) { callback(err, count); });
                            }
                        ],
                        // 'workflow_stats' is an array that has the number of 'completed','failed' and 'pending' workflows for that specfic workflow
                        function(err, workflow_stats) {
                            if (err) {
                                console.error(err);
                            } else {
			    	// Here we call the callback of the most outer layer of async.parallel, which executes 'big_task'. This callback (callback1) will
				// pass results to the callback of async.parallel whose tasks is for all workflows.
                                callback1(err,{workflow: workflow, stats: workflow_stats});
                            }
                        });
                    };
                    tasks.push(task_buffer);
	    });
            // Once we are done getting the stats for all the different workflows (i.e. 'tasks'), we will
            // get the stats for all workflows using the same method and then pass all data needed to draw this 
	    // table to the callback of async.parallel of 'big_task'
            async.parallel(tasks, 
	        function (err, each_workflow_stats) {
                    // This async.parallel is used to collect data for all workflows. 
		    async.parallel(
		        [function(callback) {
                            Report.count({workflow_run_type: 'completed'}, 
			        function(err, count) { callback(err, count); });
                         },
                         function(callback) {
                            Report.count({workflow_run_type: 'failed'}, 
			        function(err, count) { callback(err, count); });
                         },
                         function(callback) {
                            Report.count({workflow_run_type: 'pending'}, 
			        function(err, count) { callback(err, count); });
                        }],
                        function(err, all_workflow_stats) {
                            //pass the results as a two-element array to the callback of async.parallel that performs 'big_tasks'
                            callback2(null, [each_workflow_stats, all_workflow_stats]);
                        });
            });
        });
    };

    big_task.push(sum_task);

    //run four tasks in parallel, which are complete, failed, pending workflow-short_code tables and workflow-workflow_run_type table.
    async.parallel(big_task,
        function(err, results) {
            //send data to client side
            res.render('stats', {
                complete_table: results[0][0],
                complete_samples: results[0][1],
                failed_table: results[1][0],
                failed_samples: results[1][1],
                pending_table: results[2][0],
                pending_samples: results[2][1],
                stats: results[3][0],
                all: results[3][1]
            });
        });

});
module.exports = router;
