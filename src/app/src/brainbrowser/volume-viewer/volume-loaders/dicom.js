/*
* BrainBrowser: Web-based Neurological Visualization Tools
* (https://brainbrowser.cbrain.mcgill.ca)
*
* Copyright (C) 2011-2014
* The Royal Institution for the Advancement of Learning
* McGill University
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU Affero General Public License as
* published by the Free Software Foundation, either version 3 of the
* License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU Affero General Public License for more details.
*
* You should have received a copy of the GNU Affero General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
* Author: Ahmed Harmouche <ahmedharmouche92@gmail.com>
*
* Loads DICOM files for volume viewer.
*/
(function() {
  "use strict"
  
  var VolumeViewer = BrainBrowser.VolumeViewer;
  
  VolumeViewer.volume_loaders.dicom = function(description, callback) {
	  var error_message;
	  
	  if(description.dicom_files) {
        loadFiles(description.dicom_files, function(result_files) {
			var worker = new Worker("src/brainbrowser/volume-viewer/workers/dicom2nifti-worker.js");
			
			worker.addEventListener("message", function(e) {
			  var niiFromDicom = e.data[0];
			  BrainBrowser.filePool.saveFile("root", niiFromDicom.name, niiFromDicom.data);
			  VolumeViewer.volume_loaders["nifti1"]({ nii_raw: niiFromDicom.data }, callback);
			  worker.terminate();
			});
			
			worker.postMessage(result_files);
		});
	  } else {
		error_message = "invalid volume description.\n" +
        "Description must contain the property 'dicom_files'.";

		BrainBrowser.events.triggerEvent("error", { message: error_message });
		throw new Error(error_message);
    }
  };
  
  function loadFiles(fileInput, callback) {
	var files = fileInput.files;
	var loadCounter = 0;
	var result_files = [];
	
	for(var i = 0; i < files.length; i++) {
	  var f = files[i];
	  var reader = new FileReader();
		
	  reader.onload = (function(theFile) {
		return function(e1) {
		  var buffer = e1.target.result;
		  var bytes = new Uint8Array(buffer);
		  result_files.push({ "name": theFile.name, "data": bytes });
				
		  if(++loadCounter === files.length) {
		    callback(result_files);
		  }
	    }
	  })(f);
		
	  reader.readAsArrayBuffer(f);
	}
  }
}());