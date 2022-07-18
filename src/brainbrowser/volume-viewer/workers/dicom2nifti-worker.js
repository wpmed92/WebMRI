//worker for converting dicom to nifti using dicom2nifti.js
importScripts("../plugins/dicom2nifti.js");

self.addEventListener("message", function(event) {
  var dicom_files = event.data;
  
  var Dicom2NiftiModule = {
	files: dicom_files,
	passBack: function(result) {
      self.postMessage(result);
	},
	arguments: ["-z", "n", "-f", "%p_%t_%s", "-o", "/niiOut", "/dicomIn"],
	TOTAL_MEMORY: 268435456
  };
  
  dicom2nifti_run(Dicom2NiftiModule);
});