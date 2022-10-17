#

#   fast.tcl - GUI for FAST - FMRIB's Automated Segmentation Tool
#
#   Stephen Smith and Yongyue Zhang, FMRIB Image Analysis Group
#
#   Copyright (C) 2001-2002 University of Oxford
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

#{{{ setups

source [ file dirname [ info script ] ]/fslstart.tcl

set VARS(history) {}

set MC 5

catch { exec sh -c "${FSLDIR}/bin/fast 2>&1 | grep Version | awk '{ print \$3 }' -" } VERSION

#}}}
#{{{ fast

proc fast { w } {

    #{{{ vars and setup

global vars FSLDIR argc argv PWD HOME VERSION entries MC

toplevel $w
wm title $w "FAST - FMRIB's Automated Segmentation Tool - v$VERSION"
wm iconname $w "FAST"
wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm

frame $w.f

set entries($w,$MC) ""

#}}}

    tixLabelFrame $w.f.input -label "Input"
    set lfinput [ $w.f.input subwidget frame ]
    #{{{ input channels

set vars(channels) 1
tixControl $lfinput.channels -label "Number of input channels " \
        -variable vars(channels) -step 1 -min 1 -max $MC -selectmode immediate \
	-command "fast:updateinputs $w"

pack $lfinput.channels -in $lfinput -side top -anchor w -padx 5 -pady 5

#}}}
    #{{{ input image (if not in medx)

set i 0
while { $i < $MC } {

    if { $i == 0 && $argc > 0 && [ string length [ lindex $argv 0 ] ] > 0 } {
	set inputname [ remove_ext [ imglob -oneperimage [ lindex $argv 0 ] ] ]
	if { [ string first / $inputname ] == 0 || [ string first ~ $inputname ] == 0 } {
	    set entries($w,$i) $inputname
	} else {
	    set entries($w,$i) ${PWD}/$inputname
	}
	set entries($w,$MC) $entries($w,$i)
    }

    FSLFileEntry $lfinput.input$i \
	-variable entries($w,$i) \
	-pattern "IMAGE" \
	-directory $PWD \
	-label "Input image" \
	-labelwidth 23 \
	-title "Select the input image" \
	-width 40 \
	-filterhist VARS(history) \
	-command "fast:select $w"

    incr i 1
}

#}}}
    #{{{ input options

frame $lfinput.inopts

label $lfinput.inopts.typelabel -text "Image type"

set vars(type) 1
tixOptionMenu $lfinput.inopts.type -variable vars(type)
$lfinput.inopts.type add command 1 -label "T1-weighted"
$lfinput.inopts.type add command 2 -label "T2-weighted"
$lfinput.inopts.type add command 3 -label "Proton density"

pack $lfinput.inopts.typelabel $lfinput.inopts.type -in $lfinput.inopts -side left

#}}}

    tixLabelFrame $w.f.output -label "Output"
    set lfoutput [ $w.f.output subwidget frame ]
    #{{{ output image (if not in medx)

    FSLFileEntry $lfoutput.output \
            -variable entries($w,$MC) \
            -pattern "*" \
            -directory $PWD \
            -label "Output image(s) basename" \
	    -labelwidth 23 \
            -title "Select a basename for the output image(s)" \
            -width 40 \
            -filterhist VARS(history)

    pack $lfoutput.output -in $lfoutput -side top -padx 5 -pady 5 -anchor w

#}}}
    #{{{ output classes

set vars(classes) 3
tixControl $lfoutput.classes -label "Number of classes" \
        -variable vars(classes) -step 1 -min 2 -max 6 -selectmode immediate

pack $lfoutput.classes -in $lfoutput -side top -anchor w -padx 5 -pady 5

#}}}
    #{{{ output options

frame $lfoutput.outoptsA

label $lfoutput.outoptsA.label -text "Output images:"

frame $lfoutput.outoptsB

label $lfoutput.outoptsB.seglabel -text "Binary segmentation:    All classes in one image"
set vars(seg) 1
checkbutton $lfoutput.outoptsB.seg -variable vars(seg)

label $lfoutput.outoptsB.segalllabel -text "One image per class"
set vars(segall) 0
checkbutton $lfoutput.outoptsB.segall -variable vars(segall)

frame $lfoutput.outoptsC

label $lfoutput.outoptsC.problabel -text "Probability maps"
set vars(prob) 0
checkbutton $lfoutput.outoptsC.prob -variable vars(prob)

label $lfoutput.outoptsC.pvlabel -text "Partial volume maps"
set vars(pv) 0
checkbutton $lfoutput.outoptsC.pv -variable vars(pv)

label $lfoutput.outoptsC.restoredlabel -text "Restored input"
set vars(restored) 0
checkbutton $lfoutput.outoptsC.restored -variable vars(restored)

label $lfoutput.outoptsC.biaslabel -text "Bias field"
set vars(bias) 0
checkbutton $lfoutput.outoptsC.bias -variable vars(bias)

pack $lfoutput.outoptsA.label \
	-in $lfoutput.outoptsA -side left
pack $lfoutput.outoptsB.seglabel $lfoutput.outoptsB.seg \
	$lfoutput.outoptsB.segalllabel $lfoutput.outoptsB.segall \
	-in $lfoutput.outoptsB -side left
pack $lfoutput.outoptsC.problabel $lfoutput.outoptsC.prob \
	$lfoutput.outoptsC.pvlabel $lfoutput.outoptsC.pv \
	$lfoutput.outoptsC.restoredlabel $lfoutput.outoptsC.restored \
	$lfoutput.outoptsC.biaslabel $lfoutput.outoptsC.bias \
	-in $lfoutput.outoptsC -side left

pack $lfoutput.outoptsA $lfoutput.outoptsB $lfoutput.outoptsC -in $lfoutput -side top -anchor w -padx 5

#}}}

    #{{{ Advanced Options

collapsible frame $w.f.opts -title "Advanced options"

tixLabelFrame $w.f.opts.b.f
set advf [ $w.f.opts.b.f subwidget frame ]

#{{{ use a-priori probability map

frame $advf.apriori

set vars(apriori) 0
checkbutton $advf.apriori.yn -variable vars(apriori) -command "fast:updateapriori $w"
label $advf.apriori.label -text "Use a-priori probability maps for initialisation"

set vars(apriori_final) 0
checkbutton $advf.apriori.yn_final -variable vars(apriori_final)
label $advf.apriori.label_final -text "and for final segmentation"

pack $advf.apriori.label $advf.apriori.yn -in $advf.apriori -side left

#}}}
#{{{ initial tissue means

set vars(means) ""
FSLFileEntry $advf.means \
	-variable vars(means) \
	-pattern "*" \
	-directory $PWD \
	-label "Use file of initial tissue-type means" \
	-labelwidth 32 \
	-title "Select a file of initial tissue-type means" \
	-width 30 \
	-filterhist VARS(history)

#}}}
#{{{ 2D segmentation

frame $advf.2d

label $advf.2d.label -text "Do 2D (slice-by-slice) segmentation instead of 3D"

set vars(2d) 0
checkbutton $advf.2d.yn -variable vars(2d)

pack $advf.2d.label $advf.2d.yn -in $advf.2d -side left

#}}}

pack $advf.apriori $advf.2d -in $advf -side top -anchor w -padx 5 -pady 5

pack $w.f.opts.b.f -in $w.f.opts.b -anchor w

#}}}

    fast:updateinputs $w 0
    pack $w.f.input $w.f.output $w.f.opts -in $w.f -side top -padx 5 -pady 0 -anchor w

    #{{{ Button Frame

    frame $w.btns
    frame $w.btns.b -relief raised -borderwidth 1
 
    button $w.apply     -command "fast:apply $w" \
        -text "Go" -width 5
 
    button $w.cancel    -command "destroy $w" \
        -text "Exit" -width 5
 
    button $w.help -command "FmribWebHelp file: ${FSLDIR}/doc/fast/index.html" \
            -text "Help" -width 5

    pack $w.btns.b -side bottom -fill x
    pack $w.apply $w.cancel $w.help -in $w.btns.b -side left -expand yes -padx 3 -pady 10 -fill y
 
    pack $w.f $w.btns -expand yes -fill both -padx 5 -pady 5

#}}}
}

#}}}
#{{{ fast:updateinputs

proc fast:updateinputs { w junk } {

    global vars MC

    set lfinput [ $w.f.input subwidget frame ]
    set advf [ $w.f.opts.b.f subwidget frame ]

    set i 0
    while { $i < $MC } {
	pack forget $lfinput.input$i
	incr i 1
    }

    if { $vars(channels) == 1 } {

	pack $lfinput.inopts -in $lfinput -side top -padx 5 -pady 5 -anchor w

	$lfinput.input0.frame.label configure -text "Input image"
	pack $lfinput.input0 -in $lfinput -side top -padx 5 -pady 5 -anchor w -after $lfinput.channels

	pack $advf.means -in $advf -side top -anchor w -padx 5 -pady 5 -after $advf.apriori 

    } else {

	pack forget $lfinput.inopts

	set i [ expr $vars(channels) - 1 ]
	while { $i >= 0 } {
	    $lfinput.input$i.frame.label configure -text "Input image [ expr $i + 1 ]"
	    pack $lfinput.input$i -in $lfinput -side top -padx 5 -pady 5 -anchor w -after $lfinput.channels
	    incr i -1
	}

	pack forget $advf.means
    }
}

#}}}
#{{{ fast:updateapriori

proc fast:updateapriori { w } {

    global vars

    set vars(apriori_final) 0

    set advf [ $w.f.opts.b.f subwidget frame ]

    pack forget $advf.apriori.yn_final $advf.apriori.label_final

    if { $vars(apriori) == 1 } {
	pack $advf.apriori.label_final $advf.apriori.yn_final -in $advf.apriori -side left
    }
}

#}}}
#{{{ fast:select

proc fast:select { w dummy } {

    global vars entries MC

    set entries($w,0) [ remove_ext $entries($w,0) ]

    set entries($w,$MC) $entries($w,0)

#    if { [ string length $entries($w,$MC) ] == 0 } {
#	set entries($w,$MC) [ file rootname $entries($w,0) ]
#    }

}

#}}}
#{{{ fast:apply

proc fast:apply { w } {

    global vars entries MC

    set inlist ""
    for { set i 0 } { $i < $vars(channels) } { incr i 1 } {
	lappend inlist $entries($w,$i)
    }

    fast:proc $vars(channels) $inlist $vars(type) $entries($w,$MC) $vars(classes) $vars(seg) $vars(segall) $vars(prob) $vars(pv) $vars(restored) $vars(bias) $vars(apriori) $vars(apriori_final) $vars(means) $vars(2d)

    update idletasks
}

#}}}
#{{{ fast:proc

proc fast:proc { channels inlist type output classes seg segall prob pv restored bias apriori apriori_final means twod } {

    #{{{ setup for running fast 

global FSLDIR HOME

set output [ file rootname $output ]

#}}}
    #{{{ run command

if { $channels == 1 } {
    set thecommand "$FSLDIR/bin/fast -t$type"
} else {
    set thecommand "$FSLDIR/bin/mfast -s $channels"
}

set thecommand "$thecommand -c $classes"

if { ! $seg } {
    set thecommand "${thecommand} -n"
}

if { $segall } {
    set thecommand "${thecommand} -os"
}

if { $prob } {
    set thecommand "${thecommand} -op"
}

if { $pv } {
    set thecommand "${thecommand} -e -ov"
}

if { $restored } {
    set thecommand "${thecommand} -or"
}

if { $bias } {
    set thecommand "${thecommand} -ob"
}

if { $apriori_final } {
    set thecommand "${thecommand} -A"
} elseif { $apriori } {
    set thecommand "${thecommand} -a"
}

if { $means != "" } {
    set thecommand "${thecommand} -m $means"
}

if { $twod } {
    set thecommand "${thecommand} -2"
}

set thecommand "$thecommand -od $output $inlist"

set thecommand [ fsl:remote $thecommand ]

puts $thecommand

catch { exec sh -c $thecommand } ErrMsg

puts "$ErrMsg\nFinished"

#}}}
}

#}}}
#{{{ call GUI

wm withdraw .
fast .rename
tkwait window .rename

#}}}

