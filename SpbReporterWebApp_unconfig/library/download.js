// Importing the node modules that are used in the script
var express = require('express');
var router = express.Router();
var Report = require('../models/reports');
var http = require('http');

// Static variables that are used within the script
var url = "http://www-pde.hpc.oicr.on.ca/cgi-bin/spbreporter.dev/getReport.pl?range=week";
var workflow_run_types = ["completed","failed", "pending"];

/*
* This variable is used for the HIGHLIGHTING feature, Upon the first time the script is called, the variable have the value 'undefined'.
* However, as it goes through the script for the first time it will be set to the most recent modification time of any report within the set of reports
* thats being downloaded.  Anytime after the first time through, the reports that have their modification time be more recent than this
* variable will be tagged as being recently updated and eventually be highlighted.  Also, the most recent modification time for the batch
* of downloads everytime through will be checked with this variable, and if it is more recent than this value, this value will be reassigned
* to the most recent modification time for that download.
*/
var most_recent_modfication_time;


/*
* This script (thats called every 15 minutes) downloads the reports from the JSON script thats on the PDE server into the App's MongoDB database.
* It uses the Mongoose 'Report' model to translate each single JSON Object that's on
* the script into a 'reports' object in our database.  First it will delete all the reports in the database and then it will create
* a Report object for all the reports on the script and then save all of them into the database.
*/
var downloadReportsToDB = function(){
	var jsonArrays = [];

	// 'Get' request to the script, we will put the response from the 'get' request into the string variable called 'str'
	var request = http.request(url, function(response){
		var str = '';
		response.on('data', function(data){
			str+=data;
		});

		// The core of the script starts once the 'get' request is finished (Here)
		response.on('end', function(){
			// Sanity check to see if the data is properly downloaded from the script, this prevents the script from
			// continuing for the times when the JSON script is down
			if (isJSON(str)){
				// Deletes all reports from our DB
				Report.remove({},function(err) {
					if (err){
						console.log('There was an error deleting database entries');
					}
					else {
						var json_object = JSON.parse(str);
						jsonArrays.push(json_object.completed);
						jsonArrays.push(json_object.failed);
						jsonArrays.push(json_object.pending);
						/* At this point 'jsonArrays' is a 2 dimensional array, where at position 0 it contains an array of JSON objects
						* that represent workflow reports of the type 'completed', in position 1 its the same thing but for failed workflow reports
						* and in position 2 its the same thing but for pending workflow reports
						*/
						var tmp_date; // The variable that will be used to find the most recent mod time for batch of downloads

						// Will go through each JSON object in the 2 dimensional array and create a 'Report' object for it, and then save it into the database
						for (i=0;i<workflow_run_types.length;++i){
							var count = 0; // Used to count number of reports for the type of reports being downloaded (logging/debugging purposes)
							if (undefined !== jsonArrays[i]){
								for (size=0;size<jsonArrays[i].length;++size){
									++count;
									var tmp_json_report = jsonArrays[i][size];
									var tmp_modification_date = new Date(tmp_json_report.lmtime);

									if (tmp_date === undefined || tmp_modification_date > tmp_date){
										tmp_date = tmp_modification_date;
									}
									// This checks if the Report object should be marked as updated (Highlighted),
									// no reports will be marked as updated when the script runs for the first time
									var updated = tmp_modification_date > most_recent_modfication_time && most_recent_modfication_time!==undefined ? true : false;

									var report_progress = (tmp_json_report.hasOwnProperty('progress')) ? tmp_json_report.progress : undefined;
									// The Report object that will be saved to DB using the mongoose model that's been created (found in /models/reports.js)
									var tmp_report = new Report({
										sample_name : tmp_json_report.sample,
										workflow_name : tmp_json_report.workflow,
										workflow_version : tmp_json_report.version,
										workflow_run_id : tmp_json_report.wrun_id,
										workflow_run_type : workflow_run_types[i],
										workflow_run_status : tmp_json_report.status_cmd,
										create_time : tmp_json_report.crtime,
										last_modified_time : tmp_json_report.lmtime,
										progress : report_progress,
										create_date : new Date(tmp_json_report.crtime),
										last_modified_date : tmp_modification_date,
										recently_modified : updated
									});
									tmp_report.save(function(err){
										if (err){
											console.log('error saving to DB');
										}
									});
								}
							}
							// Log message for debugging purposes
							console.log(workflow_run_types[i]+ ": " + count);
						}
						// Here we updated the most_recent_modification_time variable if its the first time through for the script or the most recent mod time
						// for any report in this batch of downloads was more recent than the value set for the variable
						if (most_recent_modfication_time === undefined || tmp_date > most_recent_modfication_time){
							most_recent_modfication_time = tmp_date;
						}
						console.log('finished downloading new data to DB');
					}
				});
			}

		});
	response.on('error', function(){
		console.log('error downloading from script');
	});
});
	request.end();
	request.on('error', function(){
		console.log('There was an error downloading the data');
	});
}

/*
* Function that ensures the string passed in is a JSON object, used to ensure that the script above can go on without
* giving exceptions
*/
function isJSON(string){
	try {
		JSON.parse(string);
	} catch(err){
		return false;
	}
	return true;
}

// Was used when using the non dev JSON script (Script did not supply dates in millisecond format)
function stringToDateConverter(string){
	var tmpString = string.substring(0,string.lastIndexOf(".")).replace(" ", "T");
	return new Date(tmpString);
}
module.exports.downloadReportsToDB = downloadReportsToDB;
