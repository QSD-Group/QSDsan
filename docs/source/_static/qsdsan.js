/* add a copy button to code block*/
function addcopybtn(){

	// get all code elements
	var codeblocks = $( "div.highlight pre" );

	// For each element, do the following steps
	codeblocks.each(function(i) {
	
		// define a unique id for this element and add it
		var ID = "codeblock" + (i+1);
		$(this).attr('id', ID);

		// create a button that's configured for clipboard.js
		// https://clipboardjs.com/
		// point it to the text that's in this code block
		// add the button just after the text in the code block w/ jquery
		var copybtn = '<button class="btn copybtn" data-clipboard-target="#'+ID+'"><img src="_static/copybtn.png" width="15" alt="Copy to clipboard"></button>';
		   $(this).after(copybtn);

	});

	// tell clipboard.js to look for clicks that match this query
	new Clipboard('.btn');

}

$(document).ready(function () {
// Once the DOM (document object model) is loaded for the page, attach clipboard buttons
addcopybtn();
});
