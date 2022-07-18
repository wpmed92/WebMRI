# FSL.js

FSL.js is an asm.js port (using [Emscripten](https://emscripten.org)) of [BET](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/BET) (Brain Extraction Tool) and [FLIRT](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FLIRT) (FMRIB's Linear Registration Tool) of [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSL) (FMRIB Software Library).
The port is based on FSL version [3.3.11](https://fsl.fmrib.ox.ac.uk/fsldownloads/oldversions/fsl-3.3.11-sources.tar.gz)

Try it out [here](https://wpmed92.github.io/fsljs)!

## Background

Both BET and FLIRT are essential elements of any neurogimaging toolchain. The linearly registrate two volumes (FLIRT) one has to brainstrip (BET) the input volumes first to improve registration quality.
The popularity of FSL in the neuroimaging community and my curiousity in bringing computationally intensive programs to the web led me to try and port probably the two most popular tools in FSL.
Besides porting the tools I built a small volume viewer on top of [BrainBrowser](https://brainbrowser.cbrain.mcgill.ca), and integrated bet.js and flirt.js into this viewer.

## Project overview

*/color-maps:* contains the color maps used by BrainBrowser to render volumes

*/plugin-GUIs:* contains the JSONs of plugin GUIs. GUI for the plugins (bet.js, flirt.js) is automatically generated based on the JSON.

*/src:* contains the JavaScript and asm.js code of the project.

*->brainbrowser:* contains all the brainbrowser code used for this project. The main component is volume-viewer.
    
*-->/volume-viewer:* contains the rendering logic to show volumes.

*--->/plugins:* this project extends volume-viewer with a plugin system in which bet.js and flirt.js are defined. The files defined here are the outputs of the Emscripten compilation and can be thought of as the web equivalents of native command line programs.

*--->/workers:* contains the web workers needed to run the plugins on worker threads

*--->/volume-loaders:* this project extends BrainBrowser's volume-loaders with dicom.js, which is an Emscripten port of [dicom2nifti](https://github.com/icometrix/dicom2nifti). When DICOM files are requested to be loaded in the program, they are first converted to NIfTI and after that are they loaded.


## How to use

* Serve the root folder using any http server
* Navigate to index
* Open a volume by clicking on File -> Open NIfTI/Open DICOM
* To run bet.js: Tools -> Brain extractions
* To run flirt.js (need two input files): Tools -> Registration

## License

See the [LICENSE](LICENSE.md)