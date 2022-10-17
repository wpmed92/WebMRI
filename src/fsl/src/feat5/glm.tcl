#{{{ copyright and setup 

#   GLM - setup GLM design files
#
#   Stephen Smith, FMRIB Image Analysis Group
#
#   Copyright (C) 2006 University of Oxford
#
#   Part of FSL - FMRIB's Software Library
#   http://www.fmrib.ox.ac.uk/fsl
#   fsl@fmrib.ox.ac.uk
#   
#   Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
#   Imaging of the Brain), Department of Clinical Neurology, Oxford
#   University, Oxford, UK
#   
#   
#   LICENCE
#   
#   FMRIB Software Library, Release 3.3 (c) 2006, The University of
#   Oxford (the "Software")
#   
#   The Software remains the property of the University of Oxford ("the
#   University").
#   
#   The Software is distributed "AS IS" under this Licence solely for
#   non-commercial use in the hope that it will be useful, but in order
#   that the University as a charitable foundation protects its assets for
#   the benefit of its educational and research purposes, the University
#   makes clear that no condition is made or to be implied, nor is any
#   warranty given or to be implied, as to the accuracy of the Software,
#   or that it will be suitable for any particular purpose or for use
#   under any specific conditions. Furthermore, the University disclaims
#   all responsibility for the use which is made of the Software. It
#   further disclaims any liability for the outcomes arising from using
#   the Software.
#   
#   The Licensee agrees to indemnify the University and hold the
#   University harmless from and against any and all claims, damages and
#   liabilities asserted by third parties (including claims for
#   negligence) which arise directly or indirectly from the use of the
#   Software or the sale of any products based on the Software.
#   
#   No part of the Software may be reproduced, modified, transmitted or
#   transferred in any form or by any means, electronic or mechanical,
#   without the express permission of the University. The permission of
#   the University is not required if the said reproduction, modification,
#   transmission or transference is done without financial return, the
#   conditions of this Licence are imposed upon the receiver of the
#   product, and all original and amended source code is included in any
#   transmitted product. You may be held legally responsible for any
#   copyright infringement that is caused or encouraged by your failure to
#   abide by these terms and conditions.
#   
#   You are not permitted under this Licence to use this Software
#   commercially. Use for which any financial return is received shall be
#   defined as commercial use, and includes (1) integration of all or part
#   of the source code or the Software into a product for sale or license
#   by or on behalf of Licensee to third parties or (2) use of the
#   Software or any derivative of it for research with the final aim of
#   developing software products for sale or license to a third party or
#   (3) use of the Software or any derivative of it for research with the
#   final aim of developing non-software products for sale or license to a
#   third party, or (4) use of the Software to provide any service to an
#   external organisation for which payment is received. If you are
#   interested in using the Software commercially, please contact Isis
#   Innovation Limited ("Isis"), the technology transfer company of the
#   University, to negotiate a licence. Contact details are:
#   innovation@isis.ox.ac.uk quoting reference DE/1112.

source [ file dirname [ info script ] ]/fslstart.tcl

set VARS(history) {}

#}}}

#{{{ GLM GUI

proc glm { w } {
    global fmri FSLDIR USER feat_files VARS PWD HOME
 
    #{{{ main window basic setup

feat5:setupdefaults
set fmri(analysis) 2
set fmri(filtering_yn) 0
set fmri(poststats_yn) 0
set fmri(reg_yn) 0
set fmri(wizard_type) 1
set fmri(r_count) 30
set fmri(a_count) 30
set fmri(infeat) 0

toplevel $w

wm title      $w "GLM Setup"
wm iconname   $w "GLM"
wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm

tixBalloon $w.bhelp -initwait 3000 -options { background #b0ffb0 wrapLength 800 }

#}}}
    #{{{ main window widgets

tixOptionMenu $w.level -variable fmri(level)
$w.level add command 1 -label "Timeseries design"
$w.level add command 2 -label "Higher-level / non-timeseries design"

set fmri(npts) 100
tixControl $w.npts -label " # timepoints " -variable fmri(npts) \
	-step 1 -min 0 -selectmode immediate -options { entry.width 4 } -command "glm:updatenpts $w"

tixControl $w.tr -label " TR (s) " -variable fmri(tr) \
        -step 0.25 -min 0.001 -selectmode immediate -options { entry.width 4 }

tixControl $w.paradigm_hp -label " High pass filter cutoff (s) " \
        -variable fmri(paradigm_hp) -step 5 -min 1 -selectmode immediate -options { entry.width 5 }

pack $w.level $w.npts -in $w -side top -anchor w -padx 3 -pady 3

#{{{ button Frame

frame $w.btns
    
button $w.btns.wizard -command "feat5:wizard $w" -text "Wizard"

FSLFileDirSelectDialog $w.savedialog \
	-title  "Save FEAT setup" \
	-command "feat5:write $w 1 0 0"
$w.savedialog subwidget fsbox config -pattern "*.fsf" -filterhist VARS(history)
button $w.btns.save -command "$w.savedialog popup" -text "Save"

FSLFileDirSelectDialog $w.loaddialog \
	-title  "Load FEAT setup" \
	-command "source"
$w.loaddialog subwidget fsbox config -pattern "*.fsf" -filterhist VARS(history)
button $w.btns.load -command "$w.loaddialog popup" -text "Load"

button $w.btns.cancel -command "destroy $w" -text "Exit"

button $w.btns.help -command "FmribWebHelp file: ${FSLDIR}/doc/feat5/index.html" -text "Help"

pack $w.btns.wizard $w.btns.save $w.btns.load $w.btns.cancel $w.btns.help -in $w.btns -side left -expand yes

pack $w.btns -in $w -side bottom -fill x -padx 10 -pady 10

#}}}

$w.level configure -command "glm:updatelevel $w"

#}}}

    glm:updatelevel $w 0
}

#}}}
#{{{ glm:updatelevel

proc glm:updatelevel { w dummy } {
    global fmri

    pack forget $w.tr $w.paradigm_hp
    if { $fmri(level)==1 } {
	$w.npts configure -label " # timepoints "
	pack $w.tr $w.paradigm_hp -in $w -after $w.npts -side top -anchor w -padx 3 -pady 3
    } else {
	set fmri(multiple) $fmri(npts)
	$w.npts configure -label " # inputs "
    }

    if { [ info exists fmri(w_model) ] } {
	if { [ winfo exists $fmri(w_model) ] } {
	    destroy $fmri(w_model)
	}
    }

    feat5:setup_model_vars_simple $w
    feat5:setup_model $w
}

#}}}
#{{{ glm:updatenpts

proc glm:updatenpts { w dummy } {
    global fmri

    if { $fmri(level) == 2 && $fmri(multiple)!=$fmri(npts) } {
	glm:updatelevel $w 0
    }
}

#}}}

#{{{ call GUI 

if { ! [ info exists INGUI ] } {
    wm withdraw .
    glm .r
    tkwait window .r
}

#}}}

