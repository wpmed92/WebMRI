# WebMRI

WebMRI is an extensible web-based neuroimaging platform built on top of [BrainBrowser](https://brainbrowser.cbrain.mcgill.ca). It extends BrainBrowser VolumeViewer with a plugin system, which enables running processing algorithms on the loaded volumes. WebMRI comes with two tools that were ported from C++ to WebAssembly using [Emscripten](https://emscripten.org): [BET](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/BET) (Brain Extraction Tool) and [FLIRT](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FLIRT) (FMRIB's Linear Registration Tool) of [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSL) (FMRIB Software Library). The ports are based on FSL version [3.3.11](https://fsl.fmrib.ox.ac.uk/fsldownloads/oldversions/fsl-3.3.11-sources.tar.gz).

## Requirements for a plugin

1. The plugin has to be able to read input as NIfTI file(s).
2. The plugin has to generate its output as NIfTI file(s).
3. The plugin has to defining a command-line interface and parse arguments.

## Adding a plugin

1. Define the plugin in `main.config.js`:

```JavaScript
BrainBrowser.config.set("plugins", [
    ...
    {
      name: "Name of my plugin",
      title: "Title shown on the plugin dialog",
      author: "https://link/to/plugin/author/website",
      id: "myId",
      worker: "src/brainbrowser/volume-viewer/workers/my_worker.js",
      gui: "plugin-GUIs/my_menu.json"
    }
    ...
])
```

2. Add a JSON describing the GUI of the plugin, and the command line arguments it handles:

```JavaScript
[{ 
  "type" : "infile",
  "ext": "nii",
  "text" : "Input volume"
},
{
  "type": "outfile",
  "namebuild": "betted_,%0"
},
{
  "name" : "-f",
  "type" : "number",
  "min" : 0,
  "def" : 0.5,
  "max" : 1,
  "text" : "Fractional intensity thresshold"
},
{
  "name" : "-o",
  "type" : "bool",
  "text" : "Brain outline mask"
}
...
]
```
3. Add a web worker at `src/app/src/brainbrowser/volume-viewer/workers` to invoke your plugin:

```JavaScript
//In case of Emscripten ported plugins, this JavaScript is usually the wrapper script that invokes
//the WebAssembly module. See bet2.js and flirt.js for an example.
importScripts("../plugins/my_plugin.js");

self.addEventListener("message", function(event) {
  var my_files = event.data;
  
  var MyModule = {
	files: my_files,
	passBack: function(result) {
      self.postMessage(result);
	},
	arguments: ["-some", "-arguments"]
  };
  
  plugin_run(MyModule);
});
```

4. Server `WebMRI` and see your plugins in action.


## Demo

Try out WebMRI [here](https://wpmed92.github.io/WebMRI/src/app)!


## Project overview

The projects consists of `app` and `fsl`. The `app` folder contains the code of the application based on BrainBrowser. It integrates `bet2.wasm` and `flirt.wasm`.
The `fsl` folder contains the FSL 3.3.11 codebase modified so that it builds with Emscripten.

The directory structure of `app` is as follows:

*/color-maps:* It contains the color maps used by `BrainBrowser` to render volumes.

*/plugin-GUIs:* It contains the JSONs of plugin GUIs. GUIs for the plugins are automatically generated based on the JSON.

*/src:* It contains the JavaScript and asm.js code of the project.

*->brainbrowser:* It contains all the `BrainBrowser` code used for `WebMRI`. The main component is `VolumeViewer`.
    
*-->/volume-viewer:* It contains the rendering logic to show volumes.

*--->/plugins:* `WebMRI` extends `VolumeViewer` with a plugin system in which `bet2.js` and `flirt.js` are defined. The files defined here are the outputs of the Emscripten compilation and can be thought of as the web equivalents of native command line programs.

*--->/workers:* It contains the web workers needed to run the plugins on worker threads.

*--->/volume-loaders:* `WebMRI` extends BrainBrowser's volume-loaders with `dicom.js`, which is an Emscripten port of [dicom2nifti](https://github.com/icometrix/dicom2nifti). When DICOM files are requested to be loaded in the program, they are first converted to NIfTI and after that are they loaded.


## How to use

* Serve the root folder using any http server
* Navigate to index
* Open a volume by clicking on File -> Open NIfTI/Open DICOM
* To run bet2.js: Tools -> Brain extractions
* To run flirt.js (need two input files): Tools -> Registration

## License

See the [LICENSE](LICENSE.md)