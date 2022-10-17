#{{{ copyright and setup

#   MELODIC - Multivariate Exploratory Linear Optimized Decomposition into 
#             Independent Components
#
#   Christian Beckmann and Stephen Smith, FMRIB Image Analysis Group
#
#   Copyright (C) 2004 University of Oxford
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
 
set VERSION [ exec sh -c "${FSLDIR}/bin/melodic -V | awk '{ print \$3 }'" ]

#}}}
#{{{ main GUI

proc melodic { w } {

    #{{{ vars and setup

global vars FSLDIR INMEDX argc argv PWD VERSION entries fmri

toplevel $w
wm title $w "MELODIC - v$VERSION"
wm iconname $w "MELODIC"
wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm

frame $w.f

tixBalloon $w.bhelp -initwait 5000 -options { background \#b0ffb0 }

#set fmri(font)      [ tix option get font ]
#set fmri(bold_font) [ tix option get bold_font ]

#}}}

    tixNoteBook $w.nb
    $w.nb add data  -label "Data"
    $w.nb add filtering -label "Pre-stats"
    $w.nb add melodic -label "MELODIC"

    #{{{ Data

set fmri(dataf) [ $w.nb subwidget data ]

#{{{ multiple analyses

frame $w.multiple

set fmri(multiple) 1

tixControl $w.multiple.number -label "Number of analyses " \
	-variable fmri(multiple) -step 1 -min 1 -selectmode immediate \
	-options { entry.width 3 }

set fmri(level) 1
set fmri(analysis) 0
set fmri(inputtype) 2

button $w.multiple.setup -text "Select 4D data" -command "feat5:multiple_select $w 0 \"Select input data\" "

pack $w.multiple.number $w.multiple.setup -in $w.multiple -side left -padx 5

#}}}

#{{{ output directory

set fmri(outputdir) ""

FSLFileEntry $w.outputdir \
	-variable fmri(outputdir) \
	-directory $PWD \
	-label " Output directory" \
	-title "Name the output directory" \
	-width 35 \
	-filterhist VARS(history)

#}}}

frame $w.datamain
#{{{ npts & ndelete

frame $w.nptsndelete

#{{{ ndelete

set fmri(ndelete) 0

tixControl $w.ndelete -label "Delete volumes " -variable fmri(ndelete) \
	-step 1 -min 0 -selectmode immediate -options { entry.width 3 }

$w.bhelp bind $w.ndelete -msg "The number of initial FMRI volumes to delete before any further
processing. These should be the volumes that are not wanted because
steady-state imaging has not yet been reached - typically two or three
volumes."

#}}}

pack $w.ndelete -in $w.nptsndelete -side left

#}}}
#{{{ TR & highpass

frame $w.trparadigm_hp

#{{{ TR

set fmri(tr) 3.0
tixControl $w.tr -label "TR (s) " -variable fmri(tr) \
	-step 0.25 -min 0.001 -selectmode immediate -options { entry.width 4 }
$w.bhelp bind $w.tr -msg "The time (in seconds) between scanning successive FMRI volumes."

#}}}

pack $w.tr -in $w.trparadigm_hp -side left

#}}}
pack $w.nptsndelete $w.trparadigm_hp -in $w.datamain -side top -padx 5 -pady 3 -anchor w

pack $w.multiple $w.datamain -in $fmri(dataf) -anchor w -side top

#{{{ FSL logo

set graphpic [ image create photo -file ${FSLDIR}/tcl/fsl-logo-tiny.ppm ]
button $w.logo -image $graphpic -command "FmribWebHelp file: ${FSLDIR}/doc/index.html" -borderwidth 0
pack $w.logo -in $fmri(dataf) -anchor e -side bottom -padx 5 -pady 5

#}}}

#}}}
    #{{{ Pre-statistics processing

set fmri(filteringf) [ $w.nb subwidget filtering ]

set fmri(filtering_yn) 1

#{{{ slice timing correction

frame $w.st
set fmri(stf) $w.st

set fmri(st) 0

set fmri(st_file) ""
FSLFileEntry $w.st_file \
	-variable fmri(st_file) \
	-pattern "*" \
	-directory $PWD \
	-label "" \
	-title "Select a slice order/timings file" \
	-width 20 \
	-filterhist VARS(history)

tixOptionMenu $w.st.menu -label "Slice timing correction: " -variable fmri(st) -command "melodic2:updatest $w"
$w.st.menu add command 0 -label "None"
$w.st.menu add command 1 -label "Regular up (0, 1, 2 ... n-1)"
$w.st.menu add command 2 -label "Regular down (n-1, n-2 ... 0)"
$w.st.menu add command 3 -label "Use slice order file"
$w.st.menu add command 4 -label "Use slice timings file"

pack $w.st.menu -in $fmri(stf) -side top -side left

$w.bhelp bind $w.st -msg "Slice timing correction corrects each voxel's time-series for the fact
that later processing assumes that all slices were acquired exactly
half-way through the relevant volume's acquisition time (TR), whereas
in fact each slice is taken at slightly different times.

Slice timing correction works by using (Hanning-windowed) sinc
interpolation to shift each time-series by an appropriate fraction of
a TR relative to the middle of the TR period. It is necessary to know
in what order the slices were acquired and set the appropriate option
here.

If a slice order file is to be used (e.g., to setup interleaved slice
orderings), create a text file with a single number on each line,
where the first line states which slice was acquired first, the second
line states which slices was acquired second, etc. The first slice is
numbered 1 not 0."

#}}}
#{{{ motion correction

frame $w.mc
set fmri(mcf) $w.mc

global fmri(mc)
tixOptionMenu $w.mc.menu -label "Motion correction: " -variable fmri(mc)
$w.mc.menu add command 0 -label "None"
$w.mc.menu add command 1 -label "MCFLIRT"
set fmri(mc) 1

pack $w.mc.menu -in $fmri(mcf) -side top -side left

$w.bhelp bind $w.mc -msg "You will normally want to apply motion correction; this attempts to
remove the effect of subject head motion during the
experiment. MCFLIRT uses FLIRT (FMRIB's Linear Registration Tool)
tuned to the problem of FMRI motion correction, applying rigid-body
transformations."

#}}}
#{{{ BET brain extraction

frame $w.bet

label $w.bet.label -text "BET brain extraction"

set fmri(bet_yn) 1
checkbutton $w.bet.yn -variable fmri(bet_yn) -command "melodic2:updatebet $w"

$w.bhelp bind $w.bet -msg "This uses BET brain extraction to create a brain mask from the first
volume in the FMRI data. This is normally better than simple
intensity-based thresholding in order to remove unwanted 
voxels in FMRI data. Note that here, BET is setup to run 
in a quite liberal way so that there is very little danger 
of removing valid brain voxels."

set fmri(thresh) 10

tixControl $w.bet.thresh -label "  Threshold % " -variable fmri(thresh) -step 1 -min 0 -max 100 -selectmode immediate -options { entry.width 3 }

pack $w.bet.label $w.bet.yn $w.bet.thresh -in $w.bet -side left

#}}}
#{{{ spatial filtering

set fmri(smooth) 3
tixControl $w.smooth -label "Spatial smoothing FWHM (mm) " \
	-variable fmri(smooth) -step 1 -min 0 -selectmode immediate

$w.bhelp bind $w.smooth -msg "This determines the extent of the spatial smoothing, carried out on
each volume of the FMRI data set separately. This is 
intended to reduce noise without reducing interesting 
signal; this is successful as long as the underlying 
signal area is larger than the extent of the smoothing.

To turn off spatial smoothing simply set FWHM to 0."

#}}}
#{{{ temporal filtering

frame $w.temp

label $w.temp.label -text "Temporal filtering    "

set fmri(temphp_yn) 1
label $w.temp.hplabel -text "Highpass (s)"
checkbutton $w.temp.hp_yn -variable fmri(temphp_yn)

set fmri(templp_yn) 0

$w.bhelp bind $w.temp -msg " The high pass frequency cutoff point (in seconds) specifies 
the longest temporal period you will allow.

\"Highpass temporal filtering\" uses a local fit of a 
straight line (Gaussian-weighted within the line to give 
a smooth response) to remove low frequency artefacts. It 
is generally beneficial to remove global trends which 
will - if not removed - often dominate the decomposition."

set fmri(paradigm_hp) 100
tixControl $w.paradigm_hp -variable fmri(paradigm_hp) -step 5 -min 1 -selectmode immediate -options { entry.width 5 }

pack $w.temp.label $w.temp.hplabel $w.temp.hp_yn $w.paradigm_hp -in $w.temp -side top -side left

#}}}

pack $fmri(stf) $fmri(mcf) $w.bet $w.smooth $w.temp -in $fmri(filteringf) -anchor w -pady 1 -padx 5

#}}}
    #{{{ Melodic

set fmri(melodicf) [ $w.nb subwidget melodic ]

set fmri(varnorm) 1
checkbutton $w.varnorm -variable fmri(varnorm) -text "Variance-normalise timecourses"

$w.bhelp bind $w.varnorm -msg "When switched on, Melodic will rescale each time series so 
that the estimation is more influenced by the voxel-wise 
temporal dynamics and less by a voxels' mean signal. "

#{{{ output components

frame $w.dim

set fmri(dim_yn) 1
checkbutton $w.dim.yn -variable fmri(dim_yn) -text "Automatic dimensionality estimation" -command "melodic2:updatedim $w"

pack $w.dim.yn -in $w.dim -side left

$w.bhelp bind $w.dim -msg "In order to avoid overfitting, Melodic will attempt to estimate the number of 
components from the data using Bayesian estimators for the model 
order and use PCA to reduce the data prior to the IC estimation."

set fmri(dim) 0
tixControl $w.dim.n -label "Output components" -variable fmri(dim) -step 1 -min 1 -selectmode immediate -options { entry.width 4 }

#}}}

#{{{ thresholding

frame $w.thresh

set fmri(thresh_yn) 1
checkbutton $w.thresh.yn -variable fmri(thresh_yn) -text "Threshold IC maps" -command "melodic2:updatethresh $w"

$w.bhelp bind $w.thresh -msg "Melodic uses a mixture model approach to assign significance to individual 
voxels within a spatial map. The mixture model of a single Gaussian 
distribution (for the background noise within the spatial maps) and 
2 Gamma distributions (which model the 'active' voxels contained in 
the tails of the Gaussian) is fitted to the intensity histogram of 
the Z-transformed IC maps using a restricted EM algorithm. 

From this mixture model fit, Melodic calculates voxel-wise probabilities 
of 'activation' (as the ratio of a voxels' intensity being in the 
non-background class relative to probability of the intensity being 
background noise).
Voxels above a certain threshold level are overlayed on top of 
an example volume. The default level of 0.5 will report any voxel 
where the probability of belonging to the non-background mixtures 
exceeds the probability of the voxel belonging to the background 
noise Gaussian."

set fmri(mmthresh) "0.5"
entry $w.thresh.n -textvariable fmri(mmthresh) -width 10

pack $w.thresh.yn $w.thresh.n -in $w.thresh -side left

#}}}

set fmri(ostats) 0
checkbutton $w.ostats -variable fmri(ostats) -text "Output full stats folder"

$w.bhelp bind $w.ostats -msg "
When switched on, Melodic will save the thresholded IC 
maps and the probability maps inside a folder \/stats. 
This will substantially increase the amount of space used, 
so only switch this on if you intend to use these maps."


set fmri(report) 1
checkbutton $w.report -variable fmri(report) -text "Create web page report"


pack $w.varnorm $w.dim  $w.thresh $w.ostats $w.report -in $fmri(melodicf) -anchor w -side top

#}}}

    pack $w.nb -in $w.f -side top -padx 5 -pady 5 -anchor w

    #{{{ Button Frame

    frame $w.btns
    frame $w.btns.b -relief raised -borderwidth 1
 
    button $w.apply     -command "melodic2:apply $w keep" \
        -text "Go" -width 5
    bind $w.apply <Return> {
        [winfo toplevel %W].apply invoke}
 
    button $w.cancel    -command "melodic2:destroy $w" \
        -text "Exit" -width 5
    bind $w.cancel <Return> {
        [winfo toplevel %W].cancel invoke}
 
    button $w.help -command "FmribWebHelp file: ${FSLDIR}/doc/melodic2/index.html" \
            -text "Help" -width 5
    bind $w.help <Return> {
        [winfo toplevel %W].help invoke}

    pack $w.btns.b -side bottom -fill x
    pack $w.apply $w.cancel $w.help -in $w.btns.b \
        -side left -expand yes -padx 3 -pady 10 -fill y
 
    pack $w.f $w.btns -expand yes -fill both

#}}}
}

#}}}
#{{{ melodic2:updatest

proc melodic2:updatest { w dummy } {
    global fmri

    if { $fmri(st) < 3 } {
	pack forget $w.st_file
    } else {
	pack $w.st_file -in $fmri(stf) -side left -padx 5
    }
}

#}}}
#{{{ melodic2:updatebet

proc melodic2:updatebet { w } {

    global fmri

    if { $fmri(bet_yn) == 1 } {
	pack forget $w.bet.thresh
    } else {
	pack $w.bet.thresh -in $w.bet -side left -padx 5 -after $w.bet.yn
    }
}

#}}}
#{{{ melodic2:updatedim

proc melodic2:updatedim { w } {

    global fmri

    if { $fmri(dim_yn) == 1 } {
	pack forget $w.dim.n
    } else {
	pack $w.dim.n -in $w.dim -side left -padx 5 -after $w.dim.yn
    }
}

#}}}
#{{{ melodic2:updatethresh

proc melodic2:updatethresh { w } {

    global fmri

    if { $fmri(thresh_yn) == 0 } {
	set fmri(mmthresh) 0
	pack forget $w.thresh.n
    } else {
	set fmri(mmthresh) "0.5"
	pack $w.thresh.n -in $w.thresh -side left -padx 5 -after $w.thresh.yn
    }
}

#}}}
#{{{ melodic2:apply

proc melodic2:apply { w dialog } {
    global fmri HOME FSLDIR INMEDX entries feat_files logout
    
    for { set session 1 } { $session <= $fmri(multiple) } { incr session 1 } {

	#{{{ setup filenames

set INPUT [ remove_ext $feat_files($session) ]

set OUTPUT [ new_filename ${INPUT}.ica ]

if { ! [ file writable [ file dirname $OUTPUT ] ] } {
    set OUTPUT [ new_filename ${HOME}/[ file tail $OUTPUT ] ]
}

fsl:exec "mkdir $OUTPUT"

cd $OUTPUT

set logout ${OUTPUT}/report.log

#}}}
	#{{{ prestats

#{{{ delete volumes

set total_volumes [ exec sh -c "${FSLDIR}/bin/avwnvols $INPUT 2> /dev/null" ]

if { $fmri(ndelete) > 0 } {
    set total_volumes [ expr $total_volumes - $fmri(ndelete) ]
    fsl:exec "${FSLDIR}/bin/avwroi $INPUT prefiltered_func_data $fmri(ndelete) $total_volumes"
    set INPUT prefiltered_func_data
}

set target_vol_number [ expr $total_volumes / 2 ]
fsl:exec "${FSLDIR}/bin/avwroi $INPUT example_func $target_vol_number 1"

#}}}
#{{{ slice timing correction

if { $fmri(st) > 0 } {

    set st_opts ""
    
    switch $fmri(st) {
	2 {
	    set st_opts "--down"
	}
	3 {
	    set st_opts "--ocustom=$fmri(st_file)"
	}
	4 {
	    set st_opts "--tcustom=$fmri(st_file)"
	}
    }

    fsl:exec "${FSLDIR}/bin/slicetimer -i $INPUT --out=prefiltered_func_data_st -r $fmri(tr) $st_opts"
    set INPUT prefiltered_func_data_st
}

#}}}
#{{{ motion correction

if { $fmri(mc) != 0 } {
    fsl:exec "${FSLDIR}/bin/mcflirt -in $INPUT -out prefiltered_func_data_mcf -refvol $target_vol_number" -r
    set INPUT prefiltered_func_data_mcf
}

#}}}
#{{{ BET

if { $fmri(bet_yn) } {
    fsl:exec "${FSLDIR}/bin/betfunc $INPUT prefiltered_func_data_bet"
    set INPUT prefiltered_func_data_bet
}

#}}}
#{{{ filtering

set thecommand "${FSLDIR}/bin/ip $INPUT filtered_func_data $fmri(thresh)"

if { $fmri(smooth) > 0.01 } {
    set thecommand "$thecommand -s [ expr $fmri(smooth) / 2.355 ]"
}

if { $fmri(temphp_yn) || $fmri(templp_yn) } {

    set hp_sigma_vol -1
    if { $fmri(temphp_yn) } {
	set hp_sigma_sec [ expr $fmri(paradigm_hp) / 2.0 ]
	set hp_sigma_vol [ expr $hp_sigma_sec / $fmri(tr) ]
    }

    set lp_sigma_vol -1
    if { $fmri(templp_yn) } {
	set lp_sigma_sec 2.8
	set lp_sigma_vol [ expr $lp_sigma_sec / $fmri(tr) ]
    }

    set thecommand "$thecommand -t $hp_sigma_vol $lp_sigma_vol"

}

fsl:exec "$thecommand"

set INPUT filtered_func_data

#}}}

fsl:exec "rm -rf prefiltered_func_data*"

#}}}
	#{{{ MELODIC

set thecommand "${FSLDIR}/bin/melodic -i $INPUT -o $OUTPUT -v --nobet --bgthreshold=$fmri(thresh) --tr=$fmri(tr)"

if { $fmri(dim_yn) == 1 } {
    set thecommand "$thecommand -d 0"
} else {
    set thecommand "$thecommand -d $fmri(dim)"
}

if { $fmri(thresh_yn) == 0 } {
    set thecommand "$thecommand --no_mm"
} else {
    set thecommand "$thecommand --mmthresh=\"$fmri(mmthresh)\""
}

if { $fmri(varnorm) == 0 } {
    set thecommand "$thecommand --vn"
}

if { $fmri(report) == 1 } {
    set thecommand "$thecommand --report"
}

if { $fmri(ostats) == 1 } {
    set thecommand "$thecommand --Ostats"
}

puts $thecommand
fsl:exec $thecommand 

#}}}

	puts "Done"

	if { [ file exists report/00index.html ] } {
	    puts "To view the output report point your web browser at\n${OUTPUT}/report/00index.html\n\n"
	}
    }


    update idletasks

    if {$dialog == "destroy"} {
        melodic2:destroy $w
    }
}

#}}}
#{{{ melodic2:destroy

# Summary:      Destroys melodic dialog box
proc melodic2:destroy { w } {
    destroy $w
}

#}}}
#{{{ call GUI and wait

wm withdraw .
melodic .rename
tkwait window .rename

#}}}

