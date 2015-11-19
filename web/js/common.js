/* common.js
By : MJD
Date : 2015-10-21*/
var availableTags = [];
var previousId = '';
var rectSelection = false;

function drawBlueRect(id) {
	
	// blue rect
	$("#blueRect").remove();
	
	var height = $('svg').height();
	var tagname = $("#"+id).get(0).tagName;
	var title = $("#"+id).attr("title");
	var svg = $('#container_svg').svg();	
	var svgns = "http://www.w3.org/2000/svg";
	
	// draw blue rect
	if (tagname == 'line') {
	  var width = ($("#"+id).attr("x2") - $("#"+id).attr("x1"));
	  var posX = $("#"+id).attr("x1");

	  var rect = document.createElementNS(svgns, 'rect');
        rect.setAttributeNS(null, 'x', posX);
        rect.setAttributeNS(null, 'y', 0);
        rect.setAttributeNS(null, 'height', height);
        rect.setAttributeNS(null, 'width', width);
        rect.setAttributeNS(null, 'fill', 'blue');
		rect.setAttributeNS(null, 'style', 'opacity:0.3');
		rect.setAttributeNS(null, 'id', 'blueRect');
		rect.setAttributeNS(null, 'title', title);
		
	} else { //tagname == 'rect'
	  var width = $("#"+id).attr("width");
	  var posX = $("#"+id).attr("x");

	  var rect = document.createElementNS(svgns, 'rect');
        rect.setAttributeNS(null, 'x', posX);
        rect.setAttributeNS(null, 'y', 0);
        rect.setAttributeNS(null, 'height', height);
        rect.setAttributeNS(null, 'width', width);
        rect.setAttributeNS(null, 'fill', 'blue');
		rect.setAttributeNS(null, 'style', 'opacity:0.3');
		rect.setAttributeNS(null, 'id', 'blueRect');
		rect.setAttributeNS(null, 'title', title);
	}
   document.getElementById("group").appendChild(rect);
	
}

// zoom in image...
function zoomIn(newScale) {
	var width = $("#zoom_width_value").val();
	var height = $("#zoom_height_value").val();
	var position_x = $("#position_x").val();

	if (newScale <= 10) {
		
		var scale = $("#zoom_scale_value").val();
		var newWidth = parseInt(width) * newScale;
		var newHeight = parseInt(height) * newScale;

		// save newScale + newWidth
		$("#zoom_scale_value").val(newScale);
		
		//change scale on svg....
		var g = document.getElementById("group");
		g.setAttribute('transform','scale(' + newScale + ')');
		
		$("svg").attr("width",newWidth);
		$("svg").attr("height",newHeight);
			
		//animate "scroll-content" to "ID"
		var new_position_x = (position_x / scale) * newScale;
		$( ".scroll-content" ).animate({ marginLeft: -new_position_x}, 0);
		
		var remainder = $( ".scroll-pane" ).width() - $( ".scroll-content" ).width();
		var leftVal = parseInt( -new_position_x ); 
		var percentage = Math.round( leftVal / remainder * 100 );
		$( ".scroll-bar" ).slider( "value", percentage );
		
		$("#position_x").val( new_position_x );
			
		// zoom in
		if (newScale > 1) {
			$("#zoom-out").removeClass("disabled");
		}
		
		if (newScale >= 10) {
			$("#zoom-in").addClass("disabled");
		}
	}
}

// zoom out image ...
function zoomOut(scale) {
	var width = $("#zoom_width_value").val();
	var height = $("#zoom_height_value").val();
	var position_x = $("#position_x").val();
		
	if (scale > 1) {
		var newScale = parseInt(scale) - 1;
		var newWidth = parseInt(width) * newScale;
		var newHeight = parseInt(height) * newScale;
			
		// save newScale + newWidth
		$("#zoom_scale_value").val(newScale);

		//change scale on svg....
		var g = document.getElementById("group");
		g.setAttribute('transform','scale(' + newScale + ')');
		
		$("svg").attr("width",newWidth);
		$("svg").attr("height",newHeight);
			
		//animate "scroll-content" to "ID"
		var new_position_x = (position_x / scale) * newScale;

		var moveXLimit = ($( ".scroll-content" ).width() - $(".scroll-pane").width() + 20);

		if ( new_position_x <= 0 ) {
			new_position_x = 0; // do not center
		} else if (new_position_x > moveXLimit) {
			new_position_x = moveXLimit; // do not center
		}
		
		$( ".scroll-content" ).animate({ marginLeft: -new_position_x}, 0);
		
		var remainder = $( ".scroll-pane" ).width() - $( ".scroll-content" ).width();
		var leftVal = parseInt( -new_position_x ); 
		var percentage = Math.round( leftVal / remainder * 100 );
		$( ".scroll-bar" ).slider( "value", percentage );
		
		$("#position_x").val( new_position_x );
			
		// zoom in
		if (newScale == 1) {
			$("#zoom-out").addClass("disabled");
		}
		
		if (newScale < 10) {
			$("#zoom-in").removeClass("disabled");
		}
	}
}

// refresh display svg image on mousemove
function refresh(initialized, gap, pageX) {
	// refresh only if margin-left is between 0 and - ($( ".scroll-content" ).width)
	var moveX = ( initialized.x + pageX ) - gap;
	var moveXLimit = -($( ".scroll-content" ).width() - $(".scroll-pane").width() + 20 );
	if ( (moveX <= 0) && (moveX >= moveXLimit ) ){
	
		var node = $( ".scroll-content" );

		// image move
		node.css("margin-left", moveX);
					
		// scrollbar move
		var new_position_x = parseInt(node.css("margin-left").replace('px',''));
		var remainder = $( ".scroll-pane" ).width() - $( ".scroll-content" ).width();
		var percentage = Math.round( new_position_x / remainder * 100 );
		$( ".scroll-bar" ).slider( "value", percentage );
	
		var position_x = $( ".scroll-content" ).css( "margin-left");
		position_x = position_x.replace('px', '');
		$("#position_x").val(  Math.abs(position_x) );
	}
}

//from jquery https://jqueryui.com/slider/#side-scroll
 
    //reset slider value based on scroll content position
    function resetValue() {

	  var width_content = $( ".scroll-content" ).width();
      var remainder = $( ".scroll-pane" ).width() - $( ".scroll-content" ).width();
      var leftVal = $( ".scroll-content" ).css( "margin-left" ) === "auto" ? 0 :
        parseInt( $( ".scroll-content" ).css( "margin-left" ) );
      var percentage = Math.round( leftVal / remainder * 100 );
      $( ".scroll-bar" ).slider( "value", percentage );
	}
 
    //if the slider is 100% and window gets larger, reveal content
    function reflowContent() {
		  
      var showing = $( ".scroll-content" ).width() + parseInt( $( ".scroll-content" ).css( "margin-left" ), 10 );
      var gap = $( ".scroll-pane" ).width() - showing;
      if ( gap > 0 ) {
        $( ".scroll-content" ).css( "margin-left", parseInt( $( ".scroll-content" ).css( "margin-left" ), 10 ) + gap );
      }
	  
    }

 
// from jquery https://jqueryui.com/slider/#side-scroll
function readSvgFile() {
	
	var option_select = $( "#select_svg" ).val();
	loadSvgFiles(option_select,"init");
	testSvgGroupId();
	
}

// load svg file in container (onChange)
function loadSvgFiles(option_select,type) {
	$("body").addClass("loading");
	
	$('#search_id').val("");
	previousId = '';
	availableTags = [];

	var svg = $('#container_svg').svg('get');
	
	if ( typeof svg != 'undefined' ) {
		$('#container_svg').svg('destroy');
	} 
	
	var svg = $('#container_svg').svg();
    svg.load(option_select, {addTo: true}); 
	
	$( ".scroll-content" ).css("margin-left",0);
	if (type != "init") {
		resetValue();
	}
	
	$("#position_x").val("0");	
	
    $("#zoom_scale_value").val("1");
	$("#zoom-out").addClass("disabled");
	$("#zoom-in").removeClass("disabled");

}

// function zoom in on load...
function zoomInLoad() {
	$("#zoom_width_value").val($('svg').width());
	$("#zoom_height_value").val($('svg').height());
				
	// adjust "svg height" to "scroll-pane height"
	var zoom = 2;
		
	while(  ( $('svg').height() * zoom ) < $( ".scroll-pane" ).height() ) {
		zoom = zoom + 1;
	}
	zoom = zoom + 1; //zoom to "max"...
	
	zoomIn(zoom);	
	
	$("body").removeClass("loading");

}

// test group id load...
function testSvgGroupId() 
{ 
  var result = document.getElementById('group'); 
  var busy = false; 
  var processor = setInterval(function() 
  { 
    if(!busy) 
    { 
      busy = true; 
	  
	  result = document.getElementById('group'); 
      if (result != null) {
	  	clearInterval(processor); 
		// Parse Tags
		parseSetTagsSvg();
		// Zoom in
		zoomInLoad();
	  }

      busy = false; 

    } 

  }, 100); 

}

// parse svg id to create "autocomplete field" + set "specials" tags on line[id]...
function parseSetTagsSvg() {

	$('svg line').each(function(index) {
		var id = $(this).attr('id');
		if (typeof id != 'undefined') {
			var a = availableTags.indexOf(id); 
			if (a == -1) { //tag is unique?
				availableTags.push(id);
			}
		}
			
	});

	$('svg rect').each(function(index) {
		var id = $(this).attr('id');
		if (typeof id != 'undefined') {
			var a = availableTags.indexOf(id); 
			if (a == -1) { //tag is unique?
				availableTags.push(id);
			}
		}
			
	});

	
	//sort the tags
	availableTags.sort();
	
	$( "#id_search" ).autocomplete({
		maxResults: 40,
		source: function(request, response) {
			var results = $.ui.autocomplete.filter(availableTags, request.term);
			response(results.slice(0, this.options.maxResults));
		},
		select: function( event, ui ) {
		  if(ui.item){
            $(event.target).val(ui.item.value);
          }
          // search by id...
		  search_byid();
		}
	});
		
	// set onmouseover && onmouseout on each element "svg" line with an id
	$('svg line[id]').attr('onmouseover', "evt.target.setAttribute('opacity', '0.85');");
	$('svg line[id]').attr('onmouseout', "evt.target.setAttribute('opacity', '1');");
	$('svg line[id]').attr('onclick', "drawBlueRect($(this).attr('id'))");

	// set onmouseover && onmouseout on each element "svg" rect with an id
	$('svg rect[id]').attr('onmouseover', "evt.target.setAttribute('opacity', '0.85');");
	$('svg rect[id]').attr('onmouseout', "evt.target.setAttribute('opacity', '1');");
	$('svg rect[id]').attr('onclick', "drawBlueRect($(this).attr('id'))");
	
}

function search_byid() {
	var svg = $('#container_svg').svg('get');
				
	var g = svg.getElementById($("#id_search").val());
	var id = $("#id_search").val();
	
	// show blue rect
	drawBlueRect(id);
				
	// obtain position of "ID"
	var scale = $("#zoom_scale_value").val();
		
	var tagname = $("#"+id).get(0).tagName;
	if ( tagname == 'line' ) {
	  var x = (g.getAttribute('x1')*1)*scale; //with line, attribute => x1
	} else { //tagname == 'rect'
	  var x = (g.getAttribute('x')*1)*scale; //with rect, attribute => x
	  var current_width = (g.getAttribute('width')*1)*scale;
	  var current_height = (g.getAttribute('height')*1)*scale;
	}
		
	// center x
	var widthScrollPane = $( ".scroll-pane" ).width();
	x = x - (widthScrollPane/2);
		
	// check limit...
	var moveXLimit = ($( ".scroll-content" ).width() - $(".scroll-pane").width() + 20);
				
	if ( x < 0 ) {
		x = 0; // do not center
	} else if (x > moveXLimit) {
		x = moveXLimit; // do not center
	}
		
	//animate "scroll-content" to "ID"
	$( ".scroll-content" ).animate({ marginLeft: -x}, 1000);
		
	// change color or "ID"
	if ( tagname == 'line' ) {
	  $("#"+id).animate({svgStroke: "black", svgStrokeWidth: 10}, 2000);		
	  $("#"+id).animate({svgStroke: "black", svgStrokeWidth: 3}, 2000);			
		  
	} else { //tagname == 'rect'
	  g.setAttribute('fill', "black");
	  $("#"+id).animate({svgStroke: "black", svgStrokeWidth: 1.5}, 2000);		
	  $("#"+id).animate({svgStroke: "black", svgStrokeWidth: 0.1}, 2000);			
	}
		
	var remainder = $( ".scroll-pane" ).width() - $( ".scroll-content" ).width();
	var leftVal = parseInt( -x ); 
	var percentage = Math.round( leftVal / remainder * 100 );
	$( ".scroll-bar" ).slider( "value", percentage );
		
	$("#position_x").val( x );
		
	//reset previous id
	if ( (previousId != '') && ( $("#id_search").val() !=  previousId) ){
		var g = svg.getElementById(previousId);
		if ( tagname == 'line' ) {
		  g.setAttribute('stroke', 'green');
		} else { //tagname == 'rect'
		  g.setAttribute('fill', 'green');
		}
	}
	previousId = $("#id_search").val();
}

$( document ).ready(function() {
   	// read "first" svg file
	readSvgFile();
	
	// on "change"
    $("#select_svg").on("change", function () {
		var option_select = $( "#select_svg" ).val();
		loadSvgFiles(option_select,"notinit");
		testSvgGroupId();
    });

	// search button click
	$("a.search").on("click", function () {
		search_byid();
	});
	
	// "on keypress - Enter" in search field
	$( "#id_search" ).keypress(function(event) {
	  if ( event.which == 13 ) { // 13 => enter key
	    search_byid();
		event.preventDefault();
	  }
	});
	
	// button zoom-in
	$("#zoom-in").on("click", function () {
		var newScale = parseInt( $("#zoom_scale_value").val() ) + 1;
		zoomIn(newScale);
	});

	// button zoom-out
	$("#zoom-out").on("click", function () {
		var scale = $("#zoom_scale_value").val();
		zoomOut(scale);
	});
	
	// from jquery... https://jqueryui.com/slider/#side-scroll
	
		//scrollpane parts
		var scrollPane = $( ".scroll-pane" ),
			scrollContent = $( ".scroll-content" );
 
		//build slider
		var scrollbar = $( ".scroll-bar" ).slider({
			animate: true,
			slide: function( event, ui ) {
				if ( scrollContent.width() > scrollPane.width() ) {
					scrollContent.css( "margin-left", Math.round(
						ui.value / 100 * ( scrollPane.width() - scrollContent.width() - 20 ) /* 20 => "vertical" scrollbar */
						) + "px" );
					var position_x = scrollContent.css( "margin-left");
					position_x = position_x.replace('px', '');
					$("#position_x").val(  Math.abs(position_x) );
				} else {
					scrollContent.css( "margin-left", 0 );
					$("#position_x").val(0);
				}
			}
		});
 
		//append icon to handle
		var handleHelper = scrollbar.find( ".ui-slider-handle" )
			.mousedown(function() {
				scrollbar.width( handleHelper.width() );
			})
			.mouseup(function() {
				scrollbar.width( handleHelper.width() );
			})
			.append( "<span class='ui-icon ui-icon-grip-dotted-vertical'></span>" )
			.wrap( "<div class='ui-handle-helper-parent'></div>" ).parent();
  
		//change handle position on window resize
		$( window ).resize(function() {
			resetValue();
			reflowContent();
		});
	
	// from jquery https://jqueryui.com/slider/#side-scroll
	
		$('#container_svg').on('contextmenu',function(){return false;}); //disable default context menu on right click
		
		// source : http://upshots.org/actionscript/jquery-basics-of-dragging
		$("#container_svg").on('mousedown', function(e){

			// left click
			if ( e.which == 1 ) {
				var node = $( ".scroll-content" );
				var position = node.offset();
			
				var gap = ($( "body" ).width() -  $( ".scroll-pane" ).width()) / 2;
			
				var initialized = {
					x : position.left - e.pageX,
				};
				var previousPageX =  e.pageX;

				var handlers = {
					mousemove : function(e){
						if (Math.abs(previousPageX - e.pageX) > 5 ) {
							refresh(initialized, gap, e.pageX);
							previousPageX = e.pageX;
						}
					},
					mouseup : function(e){
						$(this).off(handlers);   
					}
				};
				$(document).on(handlers);
			}
			
			// right click
			var scale = $("#zoom_scale_value").val();
			if ((e.which == 3) && (scale < 10) ) {
				
				rectSelection = true; // display selection frame
				
				var nodeContent = $( ".scroll-content" );
				var nodePane = $( ".scroll-pane" );
				var positionContent = nodeContent.offset();
				var positionPane = nodePane.offset();
			
				var gap = ($( "body" ).width() -  $( ".scroll-pane" ).width()) / 2;
			
				var initialized = {
					x  : e.pageX - positionContent.left,
					x1  : e.pageX - positionPane.left,
					y1  : e.pageY - positionPane.top,
				};
 
				var handlers = {
					
					// show selection rectangle
					mousemove : function(e){
						// source : http://marco-difeo.de/2011/06/13/jquery-rectangle-selection/
						if (rectSelection) {
							
							// Store current mouseposition
							x2 = e.pageX - positionPane.left;
							y2 = e.pageY - positionPane.top;

							// Prevent the selection div to get outside of your frame
							//(x2+this.offsetleft < 0) ? selection = false : ($(this).width()+this.offsetleft < x2) ? selection = false : (y2 < 0) ? selection = false : ($(this).height() < y2) ? selection = false : selection = true;;
							// If the mouse is inside your frame resize the selection div
							if (rectSelection) {
								// Calculate the div selection rectancle for positive and negative values
								var TOP = (initialized.y1 < y2) ? initialized.y1 : y2;
								var LEFT = (initialized.x1 < x2) ? initialized.x1 : x2;
								var WIDTH = (initialized.x1 < x2) ? x2 - initialized.x1 : initialized.x1 - x2;
								var HEIGHT = (initialized.y1 < y2) ? y2 - initialized.y1 : initialized.y1 - y2;
								
								// Use CSS to place your selection div
								$("#selection").css({
									position: 'absolute',
									zIndex: 5000,
									left: LEFT + 20, // 20PX, position fix
									top: TOP + 5, // 5PX, position fix
									width: WIDTH,
									height: HEIGHT
								});
								$("#selection").show();

							}
						}					
						// source http://marco-difeo.de/2011/06/13/jquery-rectangle-selection/
					},
					mouseup : function(e){
						selection = false;
						$("#selection").hide();
					
						var node = $( ".scroll-content" );
						var position = node.offset();
						y2 = e.pageY - position.top;

						var mouseX = initialized.x;
						var x1 = initialized.x1;

						// left justify (x position, first click mouse)...
						$("#position_x").val(mouseX);
						
						//zoom svg image on posX
						var heightSelection = (initialized.y1 < y2) ? y2 - initialized.y1 : initialized.y1 - y2;
						var ratio = Math.floor($( ".scroll-pane" ).height() / heightSelection);
						
						var scale = $("#zoom_scale_value").val();					
						newScale = parseInt(scale) + parseInt(ratio);
						
						if (newScale > 10) {
							newScale = 10;
						}
						zoomIn(newScale); //  zoom in
						
						// correct the position to center after zoom in...
						var halfWidth = ($( ".scroll-pane" ).width() / 2) - 20; // 20px => approx. width of scrollbar
						var new_center_x = parseInt($("#position_x").val()) - halfWidth;
						
						// check limit
						var moveXLimit = ($( ".scroll-content" ).width() - $(".scroll-pane").width() + 20);
						if ( new_center_x < 0 ) {
							new_center_x = 0; // do not center
						} else if (new_center_x > moveXLimit) {
							new_center_x = moveXLimit; // do not center
						}

						$("#position_x").val( new_center_x );
						$( ".scroll-content" ).animate({ marginLeft: -new_center_x}, 0);
						
						var remainder = $( ".scroll-pane" ).width() - $( ".scroll-content" ).width();
						var leftVal = parseInt( -new_center_x ); 
						var percentage = Math.round( leftVal / remainder * 100 );
						$( ".scroll-bar" ).slider( "value", percentage );
						
								
						$(this).off(handlers);   
					}
				};
				$(document).on(handlers);
			}
			// e.which == 3
		});
		// source : http://upshots.org/actionscript/jquery-basics-of-dragging
		
		
		// source : http://stackoverflow.com/questions/1206203/how-to-distinguish-between-left-and-right-mouse-click-with-jquery
	

});