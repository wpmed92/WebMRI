/*
* Author: Ahmed Harmouche <ahmedharmouche92@gmail.com>
*/

BrainBrowser.VolumeViewer.modules.pluginsys = function(viewer) {
  "use strict";
  
  viewer.initPluginSystem = function(pluginMenuElem) {
    var guiDownloaderTasks = [];
    var loader = BrainBrowser.loader;
    var filePool = BrainBrowser.filePool;
	
    BrainBrowser.config.get("plugins").forEach(function(plugin) {
	  guiDownloaderTasks.push(loader.loadPluginGUIFromUrl(plugin));
    });
  
    Promise.all(guiDownloaderTasks).then(function(plugins) {
	  var successPlugins = plugins.filter(function(p) { return p.status === "success"; });
	
	  for(var i = 0; i < successPlugins.length; i++) {
	    var curPlugin = successPlugins[i].plugin;
	    var controlPanel = buildControlPanel(curPlugin);
	    var controlModal = $(controlModalTemplate.replace("%{id}%", curPlugin.id)
											     .replace("%{title}%", curPlugin.title));
	    controlModal.find(".modal-body").append(controlPanel);
		controlModal.find(".modal-footer").append("<div style='float:left'><a href='" + curPlugin.author + "' target='_blank'>plugin author</a></div>");
	    $("body").append(controlModal);
		$("#" + pluginMenuElem).append('<li><a href="#" data-toggle="modal" data-target="#' + curPlugin.id + '">' + curPlugin.name + '</a></li>');
      }
    }).catch(function(error) {
	   console.log(error.message);
    });
  }
  
  var optionTypeToTemplate = {
    "bool":   "<div class='checkbox'>\
	               <label>\
				     <input type='checkbox'>%{text}%\
				   </label>\
				 </div>",
    "number":  "<div>\
	               <label>\
				     <input type='number' min='%{min}%' max='%{max}%' value='%{def}%' step='%{step}%'>%{text}%\
				   </label>\
				</div>",
    "string": "<div>\
	              <label>\
				    <input type='text'>%{text}%\
				  </label>\
				</div>",
    "infile": "<div class='form-group'>\
			    <label>%{text}%</label>\
				    <select class='pluginFileSelect'></select>\
			    </select>\
		      </div>",
    "string[]": "<div class='form-group'>\
			      <label>%{text}%</label>\
				    <select>%{options}%</select>\
		         </div>"
  }
	
  var controlModalTemplate = "<div id='%{id}%' class='modal fade' role='dialog'>\
                               <div class='modal-dialog'>\
                                 <div class='modal-content'>\
                                   <div class='modal-header'>\
                                    <button type='button' class='close' data-dismiss='modal'>&times;</button>\
                                    <h4 class='modal-title'><i class='fa fa-cogs' aria-hidden='true'></i> %{title}%</h4>\
                                   </div>\
                                 <div class='modal-body'>\
                                 </div>\
								 <div class='modal-footer'>\
                                 </div>\
                               </div>\
                             </div>";
							  
  var optionMethods = {
    "toControlPanelItem" : function() {
	  var me = this;
	  var template = optionTypeToTemplate[me.type];
	  
	  if (!template) {
	    return null;
	  }
	  
	  template = template.replace("%{text}%", me.text);
	  
	  if (me.type === "string[]") {
		var options = "";
		
	    if (me.vals) {
		  for (var i = 0; i < me.vals.length; i++) {
	        options += "<option value='" + me.vals[i] + "'>" + me.vals[i] + "</option>";
		  }
		}
		
		template = template.replace("%{options}%", options);
	  } else if (me.type === "number") {
	      template = template.replace("%{min}%",  (me.min  !== undefined) ? me.min  : "")
							 .replace("%{max}%",  (me.max  !== undefined) ? me.max  : "")
							 .replace("%{def}%",  (me.def  !== undefined) ? me.def  : "0")
							 .replace("%{step}%", (me.step !== undefined) ? me.step : "0.1");
	  } 
		
	  var DOMElem = $(template);
	  me.input = DOMElem.find("input").first();
	  
	  if (!me.input.length) {
	    me.input = DOMElem.find("select").first();
	  }
	  
	  if (me.type === "infile") {
		if (me.optional) {
		  me.checkbox = $("<input type='checkbox' />");
		  DOMElem.find("label").prepend(me.checkbox);	  
		}
		
	    me.input.on('poolrefreshed', function(e) {
		  var options = "";
		  var fileNames = BrainBrowser.filePool.getFileNames((!me.ext || me.ext === "nii") ? "root" : "misc");
		  
		  for (var i = 0; i < fileNames.length; i++) {
		    options += "<option value='" + fileNames[i] + "'>" + fileNames[i] + "</option>";
		  }
		  
		  me.input.html($(options));
		});
	  }
	  
	  return DOMElem;
    },
    "getOptionValue" : function() {
	  var me = this;
	  
	  //bool
	  if (me.type === "bool") {
	    if (me.trueis && me.falseis) {
		  return { name: me.name, val: me.input.prop("checked") ? me.trueis : me.falseis };
	    } else {
		    if (me.name) {
		      return me.input.prop("checked") ? { name: me.name  } : null;
		  } else {
		      console.log("anonymus bool option needs to have 'trueis' and 'falseis' specified");
		  }
	    }
		
	  //number
	  } else if (me.type === "number") {
		  var val = me.input.val();
		  var floatVal = parseFloat(val);
			
		  return $.isNumeric(floatVal) ? { name: me.name, val: floatVal } : null;
		  
	  //infile
	  } else if (me.type === "infile"){
		  var fileName = me.input.val();
			
		  if (me.optional) {
		    if (me.checkbox.is(":checked")) {
		      return { name: me.name, val: fileName, postVal: me.postval };
			} else {
			    return null;
			}
		  }
			
		  return { name: me.name, val: fileName, postVal: me.postval };
		  
	  //outfile
	  } else if (me.type === "outfile") {
	      return { name: me.name, val: me.namebuild };
		  
	  } else {
	      var val = me.input.val();
	      return (val) ? { name: me.name, val: val } : null;
	  }
    } 
  }
  
  function addOptionMethods(option) {
    for(var methodName in optionMethods) {
	  option[methodName] = optionMethods[methodName]; 
    }
  }
  
  function buildControlPanel(plugin) {
	var optionsJSON = plugin.gui;
    var menuGroup = $("<form id='form-control-menu' class='form-group'></form>");
    var pluginExeButton = $("<button type='button' class='btn btn-default'>Run</button>");

    for(var i = 0; i < optionsJSON.length; i++) {
	  var option = optionsJSON[i];
	  addOptionMethods(option);
	  var panel = option.toControlPanelItem();
	  
	  if (panel) {
	    menuGroup.append(panel);
	  }
    }
	
    pluginExeButton.click(function(e) {
	  var args = [];
	  var files = [];
	  
	  for (var i = 0; i < optionsJSON.length; i++) {
		var curOption = optionsJSON[i];
	    var optVal = curOption.getOptionValue();
		  
		if (optVal) {
	      if (optVal.name) {
		    args.push(optVal.name);
	      }
		
	      if (optVal.val != null) {
		    if (curOption.type === "infile") {
		      var f = BrainBrowser.filePool.getFile((!curOption.ext || curOption.ext === "nii") ? "root" : "misc", optVal.val);
		      files.push({ name: f.name, data: new Uint8Array(f.data) });
		    }
		  
		    args.push(optVal.val.toString());
	      }
		  
		  if (optVal.postVal != null) {
		    args.push(optVal.postVal.toString());
		  }
		}
      }
	  
	  finalizeArguments(args);
	  runPluginWorker(plugin.worker, plugin.id, plugin.title, { "files": files, "args": args });
    });
  
    menuGroup.append(pluginExeButton);
	
    return menuGroup;
  }
  
  function finalizeArguments(args) {
	for (var i = 0; i < args.length; i++) {
	  //an 'outfile' whose name should be built based on other arguments
	  if (args[i].split(",").length > 1) {
        var nameBuild = args[i].split(",");
	    var finalName = "out/";
	  
	    for (var j = 0; j < nameBuild.length; j++) {
	      var name = nameBuild[j];
		
		  if (name.startsWith("%")) {
		    var namePart = name.substring(1);
		    var index = parseInt(namePart);
		    name = ($.isNumeric(index)) ? args[index] : args[args.indexOf(namePart) + 1];
		    name = name.substring(0, name.lastIndexOf("."));
		  }
		
		  finalName += name;
	    }
	  
	    args[i] = finalName;
	  }
	
	  console.log(args[i]);
	}
  }  
  
  function runPluginWorker(workerUrl, id, title, pluginData) {
    var worker = new Worker(workerUrl);
			
    worker.addEventListener("message", function(e) {
	  if (!BrainBrowser.filePool.dirExists(id)) {
	    BrainBrowser.filePool.mkdir(id, title);
	  } else {
	    BrainBrowser.filePool.clearDir(id);
	  }
	  
	  for (var i = 0; i < e.data.length; i++) {
		var fileName = e.data[i].name;
		var ext = fileName.substring(fileName.lastIndexOf(".") + 1);
	    BrainBrowser.filePool.saveFile((ext === "nii") ? id : "misc", e.data[i].name, e.data[i].data);
	  }
	  
	  viewer.clearVolumes();
	  viewer.loadVolume({ type: "nifti1",
						  nii_raw: BrainBrowser.filePool.getFileAtIndex(id, 0).data
					    }, function() {});
											
	  worker.terminate();
    });
	
	worker.addEventListener("onerror", function(error) {
      console.log(error.message);
	  worker.terminate();
	});
	
	worker.postMessage(pluginData);
  }
}