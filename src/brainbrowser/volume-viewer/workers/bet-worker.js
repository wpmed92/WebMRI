//worker for brain extracting the volume
importScripts("../plugins/bet2.js");

self.addEventListener("message", function(event) {
  var FSLBetModule = {
	  /*Called from Emscripten "postRun" to pass back the result volume to main*/
	  passBack: function (betted) {
		self.postMessage(betted);
	  },
	  
	  files: event.data.files,
	  
	  arguments: event.data.args,
	  
	  outputDirectory: "out",
	  
	  TOTAL_MEMORY: 956301312
  };
  
  bet2js_run(FSLBetModule);
});