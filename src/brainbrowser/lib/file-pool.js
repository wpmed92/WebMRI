(function() {
  "use strict"
  
  var fp = BrainBrowser.filePool = {
	pool: {},
	view: {},
	
	//Filepool
	init: function(DOMElem) {
	  fp.pool["root"] = { 
	    files: []
	  };
	  
	  fp.pool["misc"] = {
	    files: []
	  };
	  
	  fp.view = $("#" + DOMElem);
	  fp.view.append("<div id='div-workspace-root'><p class='text-info'><i class='fa fa-folder-open-o' aria-hidden='true'></i> root</p>\
	                  <div class='list-group' id='div-workspace-root-body'>" + fp.rootEmpty + "</div></div>");	
	},
	
	mkdir: function(dir, title) {
	  fp.pool[dir] = {
		files: []
      };
      fp.mkdirOnView(dir, title);	  
	},
	
	clearDir: function(dir) {
	  fp.pool[dir] = { files: [] };
      fp.clearDirOnView(dir); 
	},
	
	removeDir: function(dir) {
	  delete fp.pool[dir];
	  fp.removeDirFromView(dir);
	},
	
	dirExists: function(dir) {
	  return fp.pool[dir];	
	},
	
	saveFile: function(dir, name, data) {
	  /*If a non .nii file is saved (for now it's just .mat files), 
	   *save it but don't add it to the view*/
	  if (!fp.isNii(name)) {
	    fp.pool[dir].files.push({ "name": name, "data": data });
		console.log("a non-nifti file is saved:" + name);
	  } else {
		console.log("a nifti file is saved: " + name);
	    fp.pool[dir].files.push({ "name": name, "data": data });
	    fp.appendToPoolMenu(dir, name);
	  }
	  
	  if(dir === "root" || dir === "misc") {
	    fp.onPoolRefreshed();
	  }
	},
	
	getFile: function(dir, name) {
	  for (var i = 0; i < fp.pool[dir].files.length; i++) {
	    if (fp.pool[dir].files[i].name === name) {
		  return fp.pool[dir].files[i];
		}
	  } 
	},
	
	getFileAtIndex: function(dir, index) {
	  return fp.pool[dir].files[index];
	},
	
	getFileNames: function(dir) {
	  return fp.pool[dir].files.map(function(f) { return f.name; });
	},
	
	removeFile: function(dir, name) {
	  var f = fp.getFile(dir, name);
	  var index = fp.pool[dir].files.indexOf(f);
	  fp.pool[dir].files.splice(index, 1);
	  
	  if (fp.pool[dir].files.length === 0) {
		if (dir === "root") {
	      $("#div-workspace-root-body").html(fp.rootEmpty);
		}
	  } else {
	    fp.removeFromPoolMenu(dir, name);
	  }
	  
	  if (dir === "root") {
	    fp.onPoolRefreshed();
	  }
	  
	  fp.onFileRemoved();
	},
	
	moveToRoot: function(dir, name) {
	  var f = fp.getFile(dir, name);
	  var index = fp.pool[dir].files.indexOf(f);
	  
	  if (fp.pool[dir].files.length === 1) {
	    fp.removeDir(dir);
	  } else {
		fp.pool[dir].files.splice(index, 1);
	    fp.removeFromPoolMenu(dir, name);
	  }
	  
	  fp.pool['root'].files.push(f);
	  fp.appendToPoolMenu("root", name);
	  fp.onPoolRefreshed();
	},
	
	clear: function() {
	  fp.pool = {};
	  fp.clearView();
	},
	
	isNii: function(fileName) {
	  return fileName.substring(fileName.lastIndexOf('.') + 1) === "nii";
	},
	
	//Events
	onFileClick: null,
	onFileRemoved: null,
	
	onPoolRefreshed: function() {
	  $(".pluginFileSelect").trigger("poolrefreshed");
	},
	
	//View
	appendToPoolMenu: function(dir, name) {
	  var htmlName = name.replace(/\.[^/.]+$/, "-").replace(/\+/g, "plus");
	  var menuItem = $("<a id='a-" + dir + "-" + htmlName + "'class='list-group-item'><div class='a-file-pool'>" + name + "</div></a>");
	  
	  if (dir === "root") {
	    var delItem = $("<button type='button' role='button' class='btn btn-danger btn-xs pull-right'><i class='fa fa-times' aria-hidden='true'></i></button>");
		delItem.on("click", function(event) {
		  event.stopPropagation();
		  fp.removeFile(dir, name);
		});
	    menuItem.append(delItem);
	  } else {
	      var moveItem = $("<button type='button' role='button' class='btn btn-info btn-xs pull-right'><i class='fa fa-level-up' aria-hidden='true'></i></button>");
		  moveItem.on("click", function(event) {
			event.stopPropagation();
	        fp.moveToRoot(dir, name);
	      });
	      menuItem.append(moveItem);
	  }
	  
	  menuItem.on("click", function() {
	    fp.onFileClick(fp.getFile(dir, name).data);
	  });
	  
	  if (dir === "root") {
	    $("#p-root-empty").remove();
	  }
	  
	  fp.view.find("#div-workspace-" + dir + "-body").append(menuItem);
	},
	
	removeFromPoolMenu: function(dir, name) {
	  var htmlName = name.replace(/\.[^/.]+$/, "-").replace(/\+/g, "plus");
	  var id = "#a-" + dir + "-" + htmlName;
	  $(id).remove();
	},
	
	mkdirOnView: function(dir, title) {
	  var dirHeader = "<p class='text-info'><i class='fa fa-cogs' aria-hidden='true'></i> " + title + "</p>";
	  fp.view.append("<div id='div-workspace-" + dir + "'>" + dirHeader + "<div class='list-group' id='div-workspace-" + dir + "-body'></div></div>");	
	},
	
	rootEmpty: "<p id='p-root-empty'>No volumes added to root</p>",
	
	clearDirOnView: function(dir) {
	  fp.view.find("#div-workspace-" + dir + "-body").html("");
	},
	
	removeDirFromView: function(dir) {
	  fp.view.find("#div-workspace-" + dir).remove();
	},
	
	clearView: function() {
	  fp.view.html("");
	}
  }
})();