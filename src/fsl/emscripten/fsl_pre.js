function fsl_run(opts) {
	
	var isNode = typeof(exports) !== 'undefined';
	
	if (!isNode) {
		var Module = {};
		
		for (var i in opts) {
		  Module[i] = opts[i];
		}
		
		Module['preRun'] = function() {
		  //Setting the environment variable FSLOUTPUTTYPE which is normally set by ${FSLDIR}/etc/fslconf/fsl.sh 
		  ENV.FSLOUTPUTTYPE = "NIFTI";
		  FS.mkdir('/' + Module['outputDirectory']);
		  
		  for(var i = 0; i < Module['files'].length; i++) {
		    var file = Module['files'][i];
		    FS.createDataFile("/", file.name, file.data, true, true);
		  }
		};
		
		Module['postRun'] = function() {
		  Module['passBack'](getOutputFiles());
		};
		
		function getOutputFiles() {
		  var handle = FS.lookupPath(Module['outputDirectory'])
		  var buffers = [];
		  if (handle && handle.node && handle.node.contents) {
			for (var fileName in handle.node.contents) {
			  buffers.push({
			    name: fileName,
			    data: FS.readFile(Module['outputDirectory'] + "/" + fileName, {encoding: "binary"}).buffer
			  });
			}
		  }
		  return buffers;
		}
	}