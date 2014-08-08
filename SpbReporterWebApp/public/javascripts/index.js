$(document).ready(function(){
	var hash = window.location.hash;

	//---------- Code that takes care of what list to show based on what URL that was accessed-------------------

	if (hash == '#completed' || hash == ''){
		showCompleted();
	}
	else if (hash == '#pending'){
		showPending();
	}

	else if(hash == '#failed'){
		showFailed();
	}

	//---------------------- Tab selection event handling code-------------------------------------

	$('li.completed').click(function(){
		showCompleted();
	});

	$('li.pending').click(function(){
		showPending();
	});

	$('li.failed').click(function(){
		showFailed();
	});

	//------------------ Button (nav bar icon) click event handling code------------------------------------

	$('.search_button').click(function(){
		$('form.searching').toggle();
		$('.homepage-title').toggle();
		$(this).toggleClass('white');
	});

	$('.refresh_button').click(function(){
		location.reload();
	});

	$('.stats_button').click(function(){
		window.location.href = '/stats';
	});

	//------------ Code that takes care of dynamic sorting-----------------------------------

	$('input[name="sortingRadios"]').change(function(){
		var sortBy = $(this).attr('value');
		var completedItems;
		var failedItems;
		var pendingItems;

		// If the lists aren't empty put the list items (reports) in an array
		// Only sorts lists if there are more than 1 item in the list
		if (!($('.reportsList#completedReports li').hasClass('emptyReport'))){
			completedItems = $('.reportsList#completedReports li').get();
		}

		if (!($('.reportsList#failedReports li').hasClass('emptyReport'))){
			failedItems = $('.reportsList#failedReports li').get();
		}
		if (!($('.reportsList#pendingReports li').hasClass('emptyReport'))){
			pendingItems = $('.reportsList#pendingReports li').get();
		}
		//--------------------Code that takes care of sorting by date------------------
		if (sortBy == 'date'){
			if (completedItems !== undefined && completedItems.length>1){
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
			if (failedItems!== undefined && failedItems.length>1){
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
			if (pendingItems!== undefined && pendingItems.length>1){
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
			if (completedItems !== undefined && completedItems.length>1){
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
			if (failedItems !== undefined && failedItems.length>1){
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
			if (pendingItems !== undefined && pendingItems.length>1){
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
			if (completedItems !== undefined && completedItems.length>1){
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
			if (failedItems !== undefined && failedItems.length>1){
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
			if (pendingItems !== undefined && pendingItems.length>1){
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

	$('#search_form').submit(function(event){
		event.preventDefault();
		var queryString = $('input.search-div').val().toLowerCase().trim();

		if (!($('.reportsList#completedReports li').hasClass('emptyReport'))){
			$('.reportsList#completedReports li').filter(function(index){
				var workflowName = this.getElementsByClassName('workflowName')[0].textContent.toLowerCase();
				var sampleText = this.getElementsByClassName('sampleName')[0].textContent.toLowerCase();
				var sampleName = sampleText.substring(sampleText.lastIndexOf(":")+1,sampleText.length).trim();

				if($(this).hasClass('hide')){
					if (queryString == ""){
						return true;
					}
					else if (workflowName.contains(queryString) || sampleName.contains(queryString)){
						return true;
					}
					else {
						return false;
					}
				}
				else {
					if (!workflowName.contains(queryString) && !sampleName.contains(queryString)){
						return true;
					}
					else {
						return false;
					}
				}
			}).toggleClass('hide');
		}

		if (!($('.reportsList#failedReports li').hasClass('emptyReport'))){
			$('.reportsList#failedReports li').filter(function(index){
				var workflowName = this.getElementsByClassName('workflowName')[0].textContent.toLowerCase();
				var sampleText = this.getElementsByClassName('sampleName')[0].textContent.toLowerCase();
				var sampleName = sampleText.substring(sampleText.lastIndexOf(":")+1,sampleText.length).trim();

				if($(this).hasClass('hide')){
					if (queryString == ""){
						return true;
					}
					else if (workflowName.contains(queryString) || sampleName.contains(queryString)){
						return true;
					}
					else {
						return false;
					}
				}
				else {
					if (!workflowName.contains(queryString) && !sampleName.contains(queryString)){
						return true;
					}
					else {
						return false;
					}
				}
			}).toggleClass('hide');
		}

		if (!($('.reportsList#pendingReports li').hasClass('emptyReport'))){
			$('.reportsList#pendingReports li').filter(function(index){
				var workflowName = this.getElementsByClassName('workflowName')[0].textContent.toLowerCase();
				var sampleText = this.getElementsByClassName('sampleName')[0].textContent.toLowerCase();
				var sampleName = sampleText.substring(sampleText.lastIndexOf(":")+1,sampleText.length).trim();

				if($(this).hasClass('hide')){
					if (queryString == ""){
						return true;
					}
					else if (workflowName.contains(queryString) || sampleName.contains(queryString)){
						return true;
					}
					else {
						return false;
					}
				}
				else {
					if (!workflowName.contains(queryString) && !sampleName.contains(queryString)){
						return true;
					}
					else {
						return false;
					}
				}
			}).toggleClass('hide');
		}
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



