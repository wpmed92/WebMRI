#{{{ copyright and setup

#   featquery - apply masking etc to get out stats from FEAT runs
#
#   Stephen Smith, FMRIB Image Analysis Group
#
#   Copyright (C) 2002-2006 University of Oxford
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
#{{{ featquery

proc featquery { w } {
    global fmri PXHOME FSLDIR VARS argc argv PWD feat_files query vars

    toplevel $w
    wm title      $w "FEAT Stats Query"
    wm iconname   $w "FEAT Query"
    wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm

    set fmri(multiple) 1
    set fmri(level) 1
    set fmri(analysis) 4
    set fmri(inputtype) 1

    #{{{ select FEAT directories

frame $w.select -relief raised -borderwidth 1

tixControl $w.select.number -label "Number of FEAT directories " \
	-variable fmri(multiple) -step 1 -min 1 -selectmode immediate \
	-options { entry.width 3 }

button $w.select.button -text "Select" -command "feat5:multiple_select $w -1 \"Select FEAT directories\" "

pack $w.select.number $w.select.button -in $w.select -padx 5 -pady 5 -side left

#}}}
    #{{{ setup mask or coordinates

frame $w.mask -relief raised -borderwidth 1

FSLFileEntry $w.mask.mask \
	-variable fmri(mask) \
	-pattern "IMAGE" \
	-directory $PWD \
	-label "Mask image " \
	-title "Select the mask image" \
	-width 30 \
	-filterhist VARS(history)

frame $w.mask.type
set fmri(masktype) 1
radiobutton $w.mask.type.image -text "Use mask image"   -value 1 -variable fmri(masktype) -command "featquery_update $w"
radiobutton $w.mask.type.coord -text "Use co-ordinates" -value 2 -variable fmri(masktype) -command "featquery_update $w"
pack $w.mask.type.image $w.mask.type.coord -in $w.mask.type -padx 5 -side left -expand yes

pack $w.mask.mask $w.mask.type -in $w.mask -padx 5 -pady 5 -side top


frame $w.mask.coords

set fmri(coordtype) -vox
radiobutton $w.mask.coords.vox -text "vox" -value -vox -variable fmri(coordtype)
radiobutton $w.mask.coords.mm  -text "mm"  -value -mm  -variable fmri(coordtype)

set fmri(cX) 0
set fmri(cY) 0
set fmri(cZ) 0
tixControl $w.mask.coords.cX -label "  X " -variable fmri(cX) -step 1 -selectmode immediate -options { entry.width 5 }
tixControl $w.mask.coords.cY -label "Y " -variable fmri(cY) -step 1 -selectmode immediate -options { entry.width 5 }
tixControl $w.mask.coords.cZ -label "Z " -variable fmri(cZ) -step 1 -selectmode immediate -options { entry.width 5 }
pack $w.mask.coords.vox $w.mask.coords.mm $w.mask.coords.cX $w.mask.coords.cY $w.mask.coords.cZ -in $w.mask.coords -padx 5 -side left -expand yes

#}}}
    #{{{ options

frame $w.opts -relief raised -borderwidth 1

#{{{ tsplot

# frame $w.opts.tsplot

# set fmri(tsplot_yn) 1
# checkbutton $w.opts.tsplot.yn -variable fmri(tsplot_yn) -text "Create time-series plots"

# pack $w.opts.tsplot.yn -in $w.opts.tsplot -padx 5 -side left

#}}}
#{{{ percent

frame $w.opts.percent

set fmri(showpercent) 0
checkbutton $w.opts.percent.yn -variable fmri(showpercent) -text "Convert PE/COPE values to %"

pack $w.opts.percent.yn -in $w.opts.percent -padx 5 -side left

#}}}
#{{{ mask weighting

frame $w.opts.maskweight

set fmri(maskweight_yn) 0
checkbutton $w.opts.maskweight.yn -variable fmri(maskweight_yn) -text "Do not binarise mask (allow weighting)"

pack $w.opts.maskweight.yn -in $w.opts.maskweight -padx 5 -side left

#}}}
#{{{ interp thresholding

frame $w.opts.ithresh

set fmri(ithresh_yn) 0
checkbutton $w.opts.ithresh.yn -variable fmri(ithresh_yn) -text "Change post-interpolation thresholding of mask" -command "featquery_update $w"

set fmri(ithresh) 0.5
tixControl $w.opts.ithresh.thresh -variable fmri(ithresh) -step 0.1 -min 0 -max 1 -selectmode immediate -options { entry.width 5 }

pack $w.opts.ithresh.yn $w.opts.ithresh.thresh  -in $w.opts.ithresh -padx 5 -side left

#}}}
#{{{ thresholding

frame $w.opts.thresh

set fmri(statsthresh_yn) 0
checkbutton $w.opts.thresh.yn -variable fmri(statsthresh_yn) -text "Threshold stats images as well as masking" -command "featquery_update $w"

set fmri(statsthresh) 0
tixControl $w.opts.thresh.thresh -variable fmri(statsthresh) -step 1 -selectmode immediate -options { entry.width 5 }

pack $w.opts.thresh.yn $w.opts.thresh.thresh  -in $w.opts.thresh -padx 5 -side left

#}}}

#pack $w.opts.tsplot $w.opts.percent $w.opts.maskweight $w.opts.ithresh $w.opts.thresh -in $w.opts -side top -anchor w
pack $w.opts.percent $w.opts.maskweight $w.opts.ithresh $w.opts.thresh -in $w.opts -side top -anchor w

#}}}
    #{{{ output directory name

frame $w.output -relief raised -borderwidth 1

frame $w.output.name
label $w.output.name.label -text "Featquery output directory name"
set fmri(output) "featquery"
entry $w.output.name.entry -textvariable fmri(output) -width 15
pack $w.output.name.label $w.output.name.entry -in $w.output.name -padx 5 -side left

set fmri(fqpopup) 1
checkbutton $w.output.popup -variable fmri(fqpopup) -text "Popup results in web browser"

pack $w.output.name $w.output.popup -in $w.output -padx 5 -pady 5 -side top -anchor w

#}}}
    #{{{ bottom buttons

frame $w.btns -relief raised -borderwidth 1

button $w.btns.apply  -text "Go"   -command "featquery_proc $w"
button $w.btns.cancel -text "Exit" -command "destroy $w"
button $w.btns.help   -text "Help" -command "FmribWebHelp file: ${FSLDIR}/doc/feat5/featquery.html"

pack $w.btns.apply $w.btns.cancel $w.btns.help -in $w.btns -padx 5 -pady 5 -side left -expand yes

#}}}

    pack $w.select $w.mask $w.opts $w.output $w.btns -in $w -padx 5 -pady 5 -fill x 
    featquery_update $w
}

#}}}
#{{{ featquery_update

proc featquery_update { w } {
    global fmri

    pack forget $w.mask.coords
    if { $fmri(masktype) == 2 } {
	pack $w.mask.coords -in $w.mask -padx 5 -pady 5 -side top
    }

    pack forget $w.opts.thresh.thresh
    if { $fmri(statsthresh_yn) } {
	pack $w.opts.thresh.thresh -in $w.opts.thresh -after $w.opts.thresh.yn -padx 5 -side left
    }

    pack forget $w.opts.ithresh.thresh
    if { $fmri(ithresh_yn) } {
	pack $w.opts.ithresh.thresh -in $w.opts.ithresh -after $w.opts.ithresh.yn -padx 5 -side left
    }

}

#}}}
#{{{ featquery_whichstats

proc featquery_whichstats { w } {
    global fmri feat_files PWD

    if { [ winfo exists $w.stats ] } {
	destroy $w.stats
    }
    frame $w.stats -relief raised -borderwidth 1
    pack $w.stats -in $w -after $w.select -padx 5 -pady 5 -fill x

    #{{{ setup scrollbar viewport as $v

    set w0 $w.stats

    frame $w0.f
    pack $w0.f -in $w0 -side top -anchor w
    canvas $w0.f.viewport -yscrollcommand "$w0.f.ysbar set" -xscrollcommand "$w0.xsbar set" -borderwidth 0
    scrollbar $w0.xsbar -command "$w0.f.viewport xview" -orient horizontal
    scrollbar $w0.f.ysbar -command "$w0.f.viewport yview" -orient vertical
    frame $w0.f.viewport.f
    $w0.f.viewport create window 0 0 -anchor nw -window $w0.f.viewport.f
    bind $w0.f.viewport.f <Configure> "feat5:scrollform_resize $w0 $w0.f.viewport"
    pack $w0.f.viewport -side left -fill both -expand true -in $w0.f

    set v $w0.f.viewport.f

#}}}

    cd $feat_files(1)

    set row 0

    set fmri(nstats) 0
    set fmri(statslist) ""

    foreach statstype { stats/pe stats/cope stats/varcope stats/tstat stats/fstat stats/zstat stats/zfstat thresh_zstat thresh_zfstat } {
    
	set statslist [ remove_ext [ lsort -dictionary [ imglob -oneperimage ${statstype}*.* ] ] ]
    
	if { [ llength $statslist ] > 0 } {

	    label $v.${row}_0 -text "$statstype  "
	    grid $v.${row}_0 -in $v -column 0 -row $row
	    
	    set column 1

	    foreach stats $statslist {
		set fmri(dostats.$fmri(nstats)) 0
		set i [ string trimleft $stats "abcdefghijklmnopqrstuvwxyz_/" ]
		checkbutton $v.${row}_$column -variable fmri(dostats.$fmri(nstats)) -text "$i "
		grid $v.${row}_$column -in $v -column $column -row $row
		incr column 1

		set fmri(statslist) "$fmri(statslist) $stats"
		incr fmri(nstats) 1
	    }
	    
	    incr row 1
	}
	
    }

    cd $PWD
}

#}}}
#{{{ featquery_proc

proc featquery_proc { w } {
    global fmri feat_files FSLDIR

    if { ( $fmri(mask) == "" ) ||
	 ( ! [ imtest $fmri(mask) ] && ! [ imtest $feat_files(1)/$fmri(mask) ] ) } {
	puts "\n\nYou need to set the mask image!\n\n"
	return
    }

    set thecommand "${FSLDIR}/bin/featquery $fmri(multiple)"

    for { set n 1 } { $n <= $fmri(multiple) } { incr n 1 } {
	set thecommand "$thecommand $feat_files($n)"
    }

    set nstats 0
    set statslist ""
    for { set n 0 } { $n < $fmri(nstats) } { incr n 1 } {
	if { $fmri(dostats.$n) == 1 } {
	    set statslist "$statslist [ lindex $fmri(statslist) $n ]"
	    incr nstats 1
	}
    }
    set thecommand "$thecommand $nstats $statslist $fmri(output)"

    if { $fmri(showpercent) == 1 } {
    	set thecommand "$thecommand -p"
    }

    if { $fmri(statsthresh_yn) } {
	set thecommand "$thecommand -t $fmri(statsthresh)"
    }

    if { $fmri(ithresh_yn) } {
	set thecommand "$thecommand -i $fmri(ithresh)"
    }

    set thecommand "$thecommand -s"
    #    if { $fmri(tsplot_yn) } {
    #	set thecommand "$thecommand -s"
    #    }

    if { $fmri(maskweight_yn) } {
	set thecommand "$thecommand -w"
    }

    if { $fmri(fqpopup) } {
	set thecommand "$thecommand -b"
    }

    set thecommand "$thecommand $fmri(mask)"

    if { $fmri(masktype) == 2 } {
	set thecommand "$thecommand $fmri(coordtype) $fmri(cX) $fmri(cY) $fmri(cZ)"
    }

    fsl:exec $thecommand
}

#}}}
#{{{ call GUI and wait

wm withdraw .
featquery .r
tkwait window .r

#}}}
