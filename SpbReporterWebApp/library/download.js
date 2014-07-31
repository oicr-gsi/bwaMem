var express = require('express');
var router = express.Router();
var Report = require('../models/reports');
var url = "http://www-pde.hpc.oicr.on.ca/cgi-bin/spbreporter/getReport.pl?range=week";
var http = require('http');
var workflow_run_types = ["completed","failed", "pending"];
var most_recent_modfication_time;

var downloadReportsToDB = function(){
	var jsonArrays = [];
	var request = http.request(url, function(response){
		var str = '';
		response.on('data', function(data){
			str+=data;
		});

		response.on('end', function(){
			if (isJSON(str)){
				Report.remove({},function(err) {
					if (err){
						console.log('There was an error deleting database entries');
					}
					else {
						var json_object = JSON.parse(str);
						jsonArrays.push(json_object.completed);
						jsonArrays.push(json_object.failed);
						jsonArrays.push(json_object.pending);
						var tmp_date;
						for (i=0;i<workflow_run_types.length;++i){
							var count = 0;
							if (undefined !== jsonArrays[i]){
								for (size=0;size<jsonArrays[i].length;++size){
									++count;
									var tmp_json_report = jsonArrays[i][size];
									var tmp_modification_date = stringToDateConverter(tmp_json_report.lmtime);

									if (tmp_date === undefined || tmp_modification_date > tmp_date){
										tmp_date = tmp_modification_date;
									}

									var updated = tmp_modification_date > most_recent_modfication_time && most_recent_modfication_time!==undefined ? true : false;
									var report_progress = (tmp_json_report.hasOwnProperty('progress')) ? tmp_json_report.progress : undefined;
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
										create_date : stringToDateConverter(tmp_json_report.crtime),
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
								//debug tool to see if it loads correct
							console.log(workflow_run_types[i]+ ": " + count);
						}
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

function isJSON(string){
	try {
		JSON.parse(string);
	} catch(err){
		return false;
	}
	return true;
}

function stringToDateConverter(string){
	var tmpString = string.substring(0,string.lastIndexOf(".")).replace(" ", "T");
	return new Date(tmpString);
}
module.exports.downloadReportsToDB = downloadReportsToDB;
