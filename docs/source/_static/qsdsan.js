/* copy button for code blocks */
function addcopybtn(){

	// get all code blocks
	var codeblc = $('div.highlight pre');

	// for each code block
	codeblc.each(function(i) {

		// define a unique ID
		var blcid = 'codeblc' + (i+1);
		$(this).atr('id', blcid);

		// create a button that's configured for clipboard.js
		// https://github.com/zenorocha/clipboard.js
		var copybtn = '<button class="btn copybtn" data-clipboard-target="#"'+blcid+'"><img src="copybtn.png" width="10" alt="Copy to clipboard"</button>';
		$(this)/after(copybtn);
	});

	// tell clipboard.js to look for clicks that match this query
	new Clipboard('.btn');

}

// add the copy button after DOM (document object model) is loaded
$(document).ready(fuction (){
addcopytn();
});