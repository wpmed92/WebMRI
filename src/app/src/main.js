//Program entry point
$(function() {
  "use strict"
  
  BrainBrowser.VolumeViewer.start("div-brainbrowser", function(viewer) {
	var panelSize = (window.innerWidth*0.8)/3 - 20;
    var color_map_config = BrainBrowser.config.get("color_maps")[0];
	
	BrainBrowser.filePool.init("div-filepool");
	
	BrainBrowser.filePool.onFileRemoved = function() {
	  viewer.clearVolumes();	
	}
	
	BrainBrowser.filePool.onFileClick = function(file) {
	  viewer.clearVolumes();
	  viewer.loadVolume({
		type: "nifti1",
		nii_raw: file
	  }, function() {
	  });
	}
	
	viewer.initPluginSystem("ul-plugin-menu");
	viewer.loadDefaultColorMapFromURL(color_map_config.url, color_map_config.cursor_color);
	viewer.setDefaultPanelSize(panelSize, panelSize);
	viewer.render();
	
	viewer.addEventListener("volumesloaded", function(event) {
	  $("#ex1").slider({});
	  $("#ex1").on("slide", function(slideEvt) {
	    var value = parseFloat(slideEvt.value);
	    viewer.volumes[0].blend_ratios[0] = 1 - value;
	    viewer.volumes[0].blend_ratios[1] = value;
	    viewer.redrawVolumes();
	  });
	  
	  $(".blend-div").removeClass("hide");
	});
	
	viewer.addEventListener("volumeuiloaded", function(event) {
      var volume = event.volume;
	  
	  $("#input-windowing").slider({ 
	    min: volume.getVoxelMin(), 
	    max: volume.getVoxelMax(),
	    value: [volume.getVoxelMin(), volume.getVoxelMax()]
	  });
	  
      $("#input-windowing").on("slide", function(slideEvt) {
	    var values = slideEvt.value;
	    volume.intensity_min = parseFloat(values[0]);
	    volume.intensity_max = parseFloat(values[1]);
	    viewer.redrawVolumes();
      });
	});
	
	viewer.addEventListener("sliceupdate", function(event) {
	  var world_coords;
	  var voxel_coords;
	  var volume = event.volume;
	  
      if (BrainBrowser.utils.isFunction(volume.getWorldCoords)) {
        world_coords = volume.getWorldCoords();
        $("#world-x").val(world_coords.x.toPrecision(6));
        $("#world-y").val(world_coords.y.toPrecision(6));
        $("#world-z").val(world_coords.z.toPrecision(6));
      }

      if (BrainBrowser.utils.isFunction(volume.getVoxelCoords)) {
        voxel_coords = volume.getVoxelCoords();
        $("#voxel-i").val(parseInt(voxel_coords.i, 10));
        $("#voxel-j").val(parseInt(voxel_coords.j, 10));
        $("#voxel-k").val(parseInt(voxel_coords.k, 10));
      }
	});
	
	$("#file-nifti").on("change", function(e) {
	  viewer.clearVolumes();
	  viewer.loadVolume({
		type: "nifti1",
		nii_file: this
	  }, function() {
		   $("#div-volume-controls").removeClass("hide");
		   $(".blend-div").addClass("hide");
	  });
    });
	
	$("#file-dicom").on("change", function(e) {
	  viewer.clearVolumes();
	  viewer.loadVolume({
		type: "dicom",
		dicom_files: this
	  }, function() {
		   $("#div-volume-controls").removeClass("hide");
		   $(".blend-div").addClass("hide");
	  });
    });
	
	$("#a-nifti").click(function(e) {
	  $("#file-nifti").val("");
      $("#file-nifti").click();
	});
	
	$("#a-dicom").click(function(e) {
	  $("#file-dicom").val("");
      $("#file-dicom").click();
	});
	
	$("#a-overlay").click(function(e) {
      viewer.clearVolumes();
	  viewer.loadVolumes({volumes: [
		{
		  type: 'nifti1',
		  nii_raw: BrainBrowser.filePool.getFileAtIndex("root", 0).data
		},
		{
		  type: 'nifti1',
		  nii_raw: BrainBrowser.filePool.getFileAtIndex("root", 1).data
		}
	  ], blend: { desc1: 1, desc2: 2}});
	});
	
  });
  
});