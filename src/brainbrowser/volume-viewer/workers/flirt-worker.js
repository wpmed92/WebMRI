//worker for registering the volumes
importScripts("../plugins/flirt.js");

self.addEventListener("message", function(event) {
  var FSLFlirtModule = {
	  /*Called from Emscripten "postRun" to pass back the result volume to main*/
	  passBack: function (out) {
		self.postMessage(out);
	  },
	  
	  files: event.data.files,
	  
	  arguments: event.data.args,
	  
	  outputDirectory: "out",
	  
	  TOTAL_MEMORY: 956301312
  };
  
  flirtjs_run(FSLFlirtModule);
});