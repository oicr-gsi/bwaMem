$(document).ready(function(){
	var hash = window.location.hash;

	if (hash == '#completed' || hash == ''){
		showCompleted();
	}
	else if (hash == '#pending'){
		showPending();
	}

	else if(hash == '#failed'){
		showFailed();
	}

	$('input[name="sortingRadios"]').change(function(){
		var sortBy = $(this).attr('value');
		var completedItems;
		var failedItems;
		var pendingItems;

		if ($('.reportsList#completedReports li.report')){
			completedItems = $('.reportsList#completedReports li').get();
		}

		if ($('.reportsList#failedReports li.report')){
			failedItems = $('.reportsList#failedReports li').get();
		}
		if ($('.reportsList#pendingReports li.report')){
			pendingItems = $('.reportsList#pendingReports li').get();
		}
		if (sortBy == 'date'){
			if (completedItems.length>1){
				completedItems.sort(function(a,b){
					var timeOne = a.getElementsByClassName('modified')[0].textContent;

					var timeTwo = b.getElementsByClassName('modified')[0].textContent;

					if (timeOne < timeTwo) return 1;
					if (timeOne > timeTwo) return -1;
					return 0;
				});
				var completedUList = $('.reportsList#completedReports');
				$.each(completedItems,function(i, li){
					completedUList.append(li);
				});
			}
			if (failedItems.length>1){
				failedItems.sort(function(a,b){
					var timeOne = a.getElementsByClassName('modified')[0].textContent;
					var timeTwo = b.getElementsByClassName('modified')[0].textContent;

					if (timeOne < timeTwo) return 1;
					if (timeOne > timeTwo) return -1;
					return 0;
				});
				var failedUList = $('.reportsList#failedReports');
				$.each(failedItems,function(i, li){
					failedUList.append(li);
				});
			}
			if (pendingItems.length>1){
				pendingItems.sort(function(a,b){
					var timeOne = a.getElementsByClassName('modified')[0].textContent;
					var timeTwo = b.getElementsByClassName('modified')[0].textContent;

					if (timeOne < timeTwo) return 1;
					if (timeOne > timeTwo) return -1;
					return 0;
				});
				var pendingUList = $('.reportsList#pendingReports');
				$.each(pendingItems,function(i, li){
					pendingUList.append(li);
				});
			}
		}
		else if (sortBy == 'workflow'){
			if (completedItems.length>1){
				completedItems.sort(function(a,b){
					var workflowOne = a.getElementsByClassName('workflowName')[0].textContent;

					var workflowTwo = b.getElementsByClassName('workflowName')[0].textContent;

					if (workflowOne < workflowTwo) return -1;
					if (workflowOne > workflowTwo) return 1;
					return 0;
				});
				var completedUList = $('.reportsList#completedReports');
				$.each(completedItems,function(i, li){
					completedUList.append(li);
				});
			}
			if (failedItems.length>1){
				failedItems.sort(function(a,b){
					var workflowOne = a.getElementsByClassName('workflowName')[0].textContent;
					var workflowTwo = b.getElementsByClassName('workflowName')[0].textContent;

					if (workflowOne < workflowTwo) return -1;
					if (workflowOne > workflowTwo) return 1;
					return 0;
				});
				var failedUList = $('.reportsList#failedReports');
				$.each(failedItems,function(i, li){
					failedUList.append(li);
				});
			}
			if (pendingItems.length>1){
				pendingItems.sort(function(a,b){
					var workflowOne = a.getElementsByClassName('workflowName')[0].textContent;
					var workflowTwo = b.getElementsByClassName('workflowName')[0].textContent;

					if (workflowOne < workflowTwo) return -1;
					if (workflowOne > workflowTwo) return 1;
					return 0;
				});
				var pendingUList = $('.reportsList#pendingReports');
				$.each(pendingItems,function(i, li){
					pendingUList.append(li);
				});
			}
		}
		else if (sortBy == 'sample'){
			if (completedItems.length>1){
				completedItems.sort(function(a,b){
					var sampleOne = a.getElementsByClassName('sampleName')[0].textContent;

					var sampleTwo = b.getElementsByClassName('sampleName')[0].textContent;

					if (sampleOne < sampleTwo) return -1;
					if (sampleOne > sampleTwo) return 1;
					return 0;
				});
				var completedUList = $('.reportsList#completedReports');
				$.each(completedItems,function(i, li){
					completedUList.append(li);
				});
			}
			if (failedItems.length>1){
				failedItems.sort(function(a,b){
					var sampleOne = a.getElementsByClassName('sampleName')[0].textContent;
					var sampleTwo = b.getElementsByClassName('sampleName')[0].textContent;

					if (sampleOne < sampleTwo) return -1;
					if (sampleOne > sampleTwo) return 1;
					return 0;
				});
				var failedUList = $('.reportsList#failedReports');
				$.each(failedItems,function(i, li){
					failedUList.append(li);
				});
			}
			if (pendingItems.length>1){
				pendingItems.sort(function(a,b){
					var sampleOne = a.getElementsByClassName('sampleName')[0].textContent;
					var sampleTwo = b.getElementsByClassName('sampleName')[0].textContent;

					if (sampleOne < sampleTwo) return -1;
					if (sampleOne > sampleTwo) return 1;
					return 0;
				});
				var pendingUList = $('.reportsList#pendingReports');
				$.each(pendingItems,function(i, li){
					pendingUList.append(li);
				});
			}
		}
	});

	$('li.completed').click(function(){
		showCompleted();
	});

	$('li.pending').click(function(){
		showPending();
	});

	$('li.failed').click(function(){
		showFailed();
	});

	$('.refresh_button').click(function(){
		location.reload();
	});

	$('.stats_button').click(function(){
		window.location.href = '/stats';
	});

	function showCompleted(){
		$('#pendingReports').hide();
		$('#failedReports').hide();
		$('#completedReports').show();
		$('li.completed').css('background-color','white');
		$('li.failed').css('background-color','black');
		$('li.pending').css('background-color','black');


	}

	function showPending(){
		$('#pendingReports').show();
		$('#failedReports').hide();
		$('#completedReports').hide();
		$('li.pending').css('background-color','white');
		$('li.completed').css('background-color','black');
		$('li.failed').css('background-color','black');


	}

	function showFailed(){
		$('#pendingReports').hide();
		$('#failedReports').show();
		$('#completedReports').hide();
		$('li.failed').css('background-color','white');
		$('li.pending').css('background-color','black');
		$('li.completed').css('background-color','black');
	}

});



