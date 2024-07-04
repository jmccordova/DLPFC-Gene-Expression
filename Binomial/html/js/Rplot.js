(function($) {
    $(document).ready(function() {
	
	$('#Rplot').scianimator({
	    'images': ['images/Rplot1.png'],
	    'width': 480,
	    'delay': 1000,
	    'loopMode': 'loop'
	});
	$('#Rplot').scianimator('play');
    });
})(jQuery);
