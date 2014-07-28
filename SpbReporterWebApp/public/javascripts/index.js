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



