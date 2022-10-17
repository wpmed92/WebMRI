#{{{ copyright and setup 

#   FEAT TCL functions library
#
#   Stephen Smith, FMRIB Image Analysis Group
#
#   Copyright (C) 1999-2006 University of Oxford
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

#}}}

#{{{ feat5:setupdefaults

proc feat5:setupdefaults { } {
    global fmri FSLDIR HOME

    set fmri(version) 5.63

    # load defaults (mandatory!)
    if { [ file exists ${FSLDIR}/etc/fslconf/feat.tcl ] } {
	source ${FSLDIR}/etc/fslconf/feat.tcl
    } else {
	MxPause "error: FEAT default settings file ${FSLDIR}/etc/fslconf/feat.tcl doesn't exist!"
	exit 1
    }

    # load user-specific defaults if they exist 
    if { [ file exists ${HOME}/.fslconf/feat.tcl ] } {
	source ${HOME}/.fslconf/feat.tcl
    }

    set fmri(design_help) "This is a graphical representation of the design matrix and parameter
contrasts.

The bar on the left is a representation of time, which starts at the
top and points downwards. The white marks show the position of every
10th volume in time. The red bar shows the period of the longest
temporal cycle which was passed by the highpass filtering.

The main top part shows the design matrix; time is represented on the
vertical axis and each column is a different (real) explanatory
variable (e.g., stimulus type). Both the red lines and the black-white
images represent the same thing - the variation of the waveform in
time.

Below this is shown the requested contrasts; each row is a different
contrast vector and each column refers to the weighting of the
relevant explanatory variable. Thus each row will result in a Z
statistic image.

If F-tests have been specified, these appear to the right of the
contrasts; each column is a different F-test, with the inclusion of
particular contrasts depicted by filled squares instead of empty
ones."

    set fmri(infeat) 1

    set fmri(filtering_yn) 1
    set fmri(sh_yn) 0
    set fmri(st_file) ""

    set fmri(stats_yn) 1
    set fmri(wizard_type) 1

    set fmri(poststats_yn) 1

    set fmri(threshmask) ""

    set fmri(reg_yn) 1
    set fmri(reginitial_highres_nonlinear_yn) 0
    set fmri(reghighres_nonlinear_yn) 0
    set fmri(regstandard_nonlinear_yn) 0

    set fmri(ncopeinputs) 0
    set fmri(inputtype) 1
    set fmri(multiple) 1
    set fmri(outputdir) ""
    set fmri(relative_yn) 0
    set fmri(constcol) 0
    set fmri(con_mode) orig
    set fmri(con_mode_old) orig
    set fmri(evs_orig) 1
    set fmri(evs_real) 1
    set fmri(ncon_real) 1
    set fmri(ncon_orig) 1
    set fmri(nftests_real) 0
    set fmri(nftests_orig) 0

    set fmri(feat_filename) [ exec sh -c "${FSLDIR}/bin/tmpnam /tmp/feat" ].fsf

}

#}}}
#{{{ feat5:scrollform_resize

proc feat5:scrollform_resize { w0 viewport } {

    set MAXWIDTH  950
    set MAXHEIGHT 700

    set fwidth  [ winfo width $viewport.f ]
    if { $fwidth > $MAXWIDTH } {
	set fwidth $MAXWIDTH
	pack $w0.xsbar -fill x -in $w0
    } else {
	pack forget $w0.f.xsbar
    }

    set fheight [ winfo height $viewport.f ]
    if { $fheight > $MAXHEIGHT } {
	set fheight $MAXHEIGHT
	pack $w0.f.ysbar -side right -fill y -in $w0.f
    } else {
	pack forget $w0.f.ysbar
    }

    $viewport configure -width $fwidth -height $fheight -scrollregion [ $viewport bbox all ]
}

#}}}
#{{{ feat5:multiple_select

proc feat5:multiple_select { w which_files windowtitle } {
    global FSLDIR fmri VARS PWD feat_files unwarp_files unwarp_files_mag initial_highres_files highres_files

    #{{{ setup window

    set count 0
    set w0 ".dialog[incr count]"
    while { [ winfo exists $w0 ] } {
        set w0 ".dialog[incr count]"
    }

    toplevel $w0

    wm iconname $w0 "Select"
    wm iconbitmap $w0 @${FSLDIR}/tcl/fmrib.xbm

    wm title $w0 $windowtitle

    frame $w0.f
    pack $w0.f -expand yes -fill both -in $w0 -side top

    canvas $w0.f.viewport -yscrollcommand "$w0.f.ysbar set"
    scrollbar $w0.f.ysbar -command "$w0.f.viewport yview" -orient vertical
    frame $w0.f.viewport.f
    $w0.f.viewport create window 0 0 -anchor nw -window $w0.f.viewport.f
    bind $w0.f.viewport.f <Configure> "feat5:scrollform_resize $w0 $w0.f.viewport"
    pack $w0.f.viewport -side left -fill both -expand true -in $w0.f

#}}}
    #{{{ setup file selections

set pastevar feat_files

for { set i 1 } { $i <= $fmri(multiple) } { incr i 1 } {

    if { $which_files < 1 } {

	if { ( $fmri(level) == 1 && $fmri(analysis) != 4 && $fmri(analysis) != 0 ) || $fmri(inputtype) == 2 } {
	    FSLFileEntry $w0.filename$i -label " $i:   " -directory $PWD -title "Select input data" -width 60 -filterhist VARS(history) \
		    -variable feat_files($i) -pattern "IMAGE"
	} else {
	    FSLFileEntry $w0.filename$i -label " $i:   " -directory $PWD -title "Select input data" -width 60 -filterhist VARS(history) \
		    -variable feat_files($i) -pattern "*.feat" -dirasfile "design.fsf"
	}

    } elseif { $which_files == 1 } {
	FSLFileEntry $w0.filename$i -label " $i:   " -directory $PWD -title "Select input data" -width 60 -filterhist VARS(history) \
		-variable unwarp_files($i) -pattern "IMAGE"
	set pastevar unwarp_files
    } elseif { $which_files == 2 } {
	FSLFileEntry $w0.filename$i -label " $i:   " -directory $PWD -title "Select input data" -width 60 -filterhist VARS(history) \
		-variable unwarp_files_mag($i) -pattern "IMAGE"
	set pastevar unwarp_files_mag
    } elseif { $which_files == 3 } {
	FSLFileEntry $w0.filename$i -label " $i:   " -directory $PWD -title "Select input data" -width 60 -filterhist VARS(history) \
		-variable initial_highres_files($i) -pattern "IMAGE"
	set pastevar initial_highres_files
    } elseif { $which_files == 4 } {
	FSLFileEntry $w0.filename$i -label " $i:   " -directory $PWD -title "Select input data" -width 60 -filterhist VARS(history) \
		-variable highres_files($i) -pattern "IMAGE"
	set pastevar highres_files
    }

    pack $w0.filename$i -in $w0.f.viewport.f -side top -expand yes -fill both -padx 3 -pady 3
}

#}}}
    #{{{ setup buttons

frame $w0.btns
frame $w0.btns.b -relief raised -borderwidth 1

button $w0.pastebutton -command "feat5:multiple_paste \"Input data\" 1 $fmri(multiple) $pastevar x" -text "Paste"

button $w0.cancel -command "feat5:multiple_check $w $which_files 1 1 0; destroy $w0" -text "OK"

pack $w0.btns.b -side bottom -fill x -padx 3 -pady 5
pack $w0.pastebutton $w0.cancel -in $w0.btns.b -side left -expand yes -padx 3 -pady 3 -fill y
pack $w0.btns -expand yes -fill x

#}}}
}

#}}}
#{{{ feat5:multiple_check

proc feat5:multiple_check { w which_files load updateimageinfo dummy } {
    global FSLDIR fmri feat_files unwarp_files unwarp_files_mag initial_highres_files highres_files

    if { $which_files < 0 } {
	set load 0
	featquery_whichstats $w
    }

    if { $which_files > 0 } {

	if { !$fmri(reg_yn) } {
	    return 0
	}

	if { [ exec sh -c "${FSLDIR}/bin/avwnvols [ file rootname $fmri(regstandard) ] 2> /dev/null" ] < 1 } {
	    MxPause "Warning: standard image not valid"
	    return 1
	}
    }

    if { $fmri(multiple) < 2 } {
	set nmultiple 1
    } else {
	set nmultiple $fmri(multiple)
    }

    set AllOK 1

    for { set i 1 } { $i <= $nmultiple } { incr i 1 } {

	if { $which_files == 0 } {

	    if { ! [ info exists feat_files($i) ] } {
		set feat_files($i) ""
	    }

	    if { ! [ file exists $feat_files($i) ] } {
		set AllOK 0
	    } else {
		if { $fmri(level) == 1 && $fmri(analysis) != 4 && $fmri(analysis) != 0 } {
		    if {  $i == 1 && $updateimageinfo } {
			feat5:updateimageinfo $w $i 1
		    }
		} else {

		    if { $fmri(inputtype) == 1 } {
			if { ! [ file exists $feat_files($i)/design.fsf ] } {
			    set AllOK 0
			} else {
			    if { $load && $i == 1 && $fmri(level)==1 } {
				feat5:load $w 0 $feat_files(1)/design.fsf
				if { $fmri(level)==2 && $fmri(mixed_yn)==1 } {
				    MxPause "Error: re-running just Post-stats is not possible on existing higher-level FEAT directories that were generated using FLAME 1+2.

To change thresholding, either fully re-run FLAME 1+2 analysis with different Post-stats settings, or fully re-run using FLAME 1, after which it is possible to re-run purely Post-stats."
				    set AllOK 0
				} else {
				    MxPause "Warning: have just loaded in design information from the design.fsf
in the first FEAT directory in the list."
				}
			    }
			}
		    }
		}
	    }

	    if { $fmri(level) == 2 && $AllOK } {
		feat5:updateselect $w 0
	    }

	} elseif { $which_files == 1 } {

	    if { ! [ info exists unwarp_files($i) ] } {
		set unwarp_files($i) ""
	    }

	    if { [ exec sh -c "${FSLDIR}/bin/avwnvols [ remove_ext $unwarp_files($i) ] 2> /dev/null" ] < 1 } {
		set AllOK 0
	    }

	} elseif { $which_files == 2 } {

	    if { ! [ info exists unwarp_files_mag($i) ] } {
		set unwarp_files_mag($i) ""
	    }

	    if { [ exec sh -c "${FSLDIR}/bin/avwnvols [ remove_ext $unwarp_files_mag($i) ] 2> /dev/null" ] < 1 } {
		set AllOK 0
	    }

	} elseif { $which_files == 3 } {

	    if { ! [ info exists initial_highres_files($i) ] } {
		set initial_highres_files($i) ""
	    }

	    if { [ exec sh -c "${FSLDIR}/bin/avwnvols [ remove_ext $initial_highres_files($i) ] 2> /dev/null" ] < 1 } {
		set AllOK 0
	    }

	} elseif { $which_files == 4 } {

	    if { ! [ info exists highres_files($i) ] } {
		set highres_files($i) ""
	    }

	    if { [ exec sh -c "${FSLDIR}/bin/avwnvols [ remove_ext $highres_files($i) ] 2> /dev/null" ] < 1 } {
		set AllOK 0
	    }

	}
    }
    
    if { ! $AllOK} {
	MxPause "Warning: you haven't filled in all the relevant selections with valid
filenames - this could include image data, FEAT directories,
unwarp input files and high-resolution images."
	return 1
    }

    return 0
}

#}}}
#{{{ feat5:multiple_paste

proc feat5:multiple_paste { windowtitle xdim ydim var1name var2 } {
    global FSLDIR
    upvar $var1name var1

    #{{{ setup window

set count 0
set w0 ".dialog[incr count]"
while { [ winfo exists $w0 ] } {
    set w0 ".dialog[incr count]"
}

toplevel $w0

wm iconname $w0 "Paste"
wm iconbitmap $w0 @${FSLDIR}/tcl/fmrib.xbm

wm title $w0 "$windowtitle - ${ydim}x$xdim Paste Window"

#}}}
    #{{{ setup main panel

scrollbar $w0.xsbar -command "$w0.text xview" -orient horizontal
scrollbar $w0.ysbar -command "$w0.text yview" -orient vertical

text $w0.text -xscrollcommand "$w0.xsbar set" -yscrollcommand "$w0.ysbar set" \
  -width 60 -height 20 -wrap none 

pack $w0.xsbar -side bottom -fill x
pack $w0.ysbar -side right  -fill y

pack $w0.text -side left -expand yes -fill both

#}}}
    #{{{ setup buttons

frame $w0.buttons

button $w0.buttons.clear -command "$w0.text delete 0.0 end" -text "Clear"

button $w0.buttons.cancel -command "feat5:multiple_paste_process $w0 $xdim $ydim $var1name $var2" -text "OK"

pack $w0.buttons -before $w0.xsbar -side bottom

pack $w0.buttons.clear $w0.buttons.cancel -in $w0.buttons -side left -padx 3 -pady 3

#}}}
    #{{{ fill paste window

    for { set y 1 } { $y <= $ydim } { incr y 1 } {
	for { set x 1 } { $x <= $xdim } { incr x 1 } {
	    if { $var2 == "x" } {
		$w0.text insert end "$var1(${y})"
	    } else {
		$w0.text insert end "$var1(${var2}${y}.${x})\t"
	    }
	}
	$w0.text insert end "\n"
    }

#}}}
}

proc feat5:multiple_paste_process { w0 xdim ydim var1name var2 } {
    upvar $var1name var1

    set alltext [ concat [ $w0.text get 0.0 end ] ]

    if { [ llength $alltext ] < [ expr $xdim * $ydim ] } {
	MxPause "Not enough entries!"
	return 1
    }

    set i 0
    for { set y 1 } { $y <= $ydim } { incr y 1 } {
	for { set x 1 } { $x <= $xdim } { incr x 1 } {
	    if { $var2 == "x" } {
		set var1(${y}) [ lindex $alltext $i ]
	    } else {
		set var1(${var2}${y}.${x}) [ lindex $alltext $i ]
	    }
	    incr i 1
	}
    }

    destroy $w0
}

#}}}
#{{{ feat5:apply

proc feat5:apply { w } {
    global fmri feat_files unwarp_files unwarp_files_mag initial_highres_files highres_files FSLDIR HOME

    #{{{ update variables

foreach v { tr ndelete paradigm_hp brain_thresh smooth uncorrected voxel prob_thresh z_thresh zmin zmax } {
    $w.$v update
}

#}}}
    #{{{ write model

if { $fmri(level) > 1 || $fmri(stats_yn) || $fmri(poststats_yn) } {
    set createthemodel 1
} else {
    set createthemodel 0
}

if { [ feat5:write $w $createthemodel 1 1 $fmri(feat_filename) ] } {
    return 1
}

#}}}
    #{{{ if level>1 test for existing registration

if { $fmri(level) > 1 && $fmri(analysis) != 4 } {

    set problem 0

    for { set i 1 } { $i <= $fmri(multiple) } { incr i 1 } {

	set featdirname $feat_files($i)
	if { $fmri(inputtype) == 2 } {
	    set featdirname [ file dirname [ file dirname $featdirname ] ]
	}

	if { ! [ file exists $featdirname/reg/example_func2standard.mat ] && 
	     ! [ file exists $featdirname/example_func2standard.mat ] && 
	     ! [ file exists $featdirname/design.lev ] } {
	    set problem 1
	}
    }

    if { $problem } {
	MxPause "Registration has not been run for all of the FEAT directories that you have selected for group analysis.

Please turn on and setup registration."
	return 1
    }

}

#}}}
    #{{{ create delay if required

set thedelay ""
if { $fmri(delay) > 0.003 } {
    set thedelay "sleep [ expr int( $fmri(delay) * 60.0 * 60.0 ) ] ; "
}

#}}}

    set FSFROOT [ file rootname $fmri(feat_filename) ]

    set feat_script $FSLDIR/bin/feat
    if { [ file exists /usr/local/bin/feat_fmrib ] } {
	set feat_script /usr/local/bin/feat_fmrib
    }

    catch { exec sh -c "$thedelay $feat_script $FSFROOT" & } junk

    set fmri(feat_filename) [ exec sh -c "${FSLDIR}/bin/tmpnam /tmp/feat" ].fsf

    update idletasks
}

#}}}
#{{{ feat5:updateanalysis and updatelevel

proc feat5:updatelevel { w dummy } {
    global fmri

    if { $fmri(level) == 1 } {
	$w.mode.analysis enable 0
	$w.mode.analysis enable 7
	$w.mode.analysis enable 1
	$w.mode.analysis enable 3
	$w.mode.analysis enable 2
	$w.mode.analysis enable 4
	set fmri(analysis) 7
	set fmri(reg_yn) 1
	set fmri(r_count) 30
	set fmri(a_count) 30
    } else {
	set fmri(analysis) 6
	set fmri(reg_yn) 0
	$w.mode.analysis disable 0
	$w.mode.analysis disable 7
	$w.mode.analysis disable 1
	$w.mode.analysis disable 3
	$w.mode.analysis disable 2
	$w.mode.analysis disable 4
    }

    feat5:setup_model_vars_simple $w
    set fmri(filmsetup) 0

    feat5:updateanalysis $w 0
}

proc feat5:updateanalysis { w dummy } {
    global fmri

    set fmri(filtering_yn) [ expr   $fmri(analysis)       % 2 == 1 ]
    set fmri(stats_yn)     [ expr ( $fmri(analysis) / 2 ) % 2 == 1 ]
    set fmri(poststats_yn) [ expr ( $fmri(analysis) / 4 ) % 2 == 1 ]

    #{{{ update notebook view

if { !$fmri(filtering_yn) } {
    $w.nb pageconfigure filtering -state disabled
}  else {
    $w.nb pageconfigure filtering -state normal
}

if { !$fmri(stats_yn) } {
    $w.nb pageconfigure stats -state disabled
}  else {
    $w.nb pageconfigure stats -state normal
}

if { !$fmri(poststats_yn) } {
    $w.nb pageconfigure poststats -state disabled
} else {
    $w.nb pageconfigure poststats -state normal
    if { !$fmri(stats_yn) } {
	set fmri(ndelete) 0
    }
}

if { $fmri(level)==1 } {
    $w.nb pageconfigure reg -state normal
}  else {
    $w.nb pageconfigure reg -state disabled
}

#}}}
    #{{{ update misc and data

pack forget $w.newdir_yn
if { $fmri(analysis) == 0 || $fmri(analysis) == 4 } {
    pack $w.newdir_yn -in $fmri(miscf) -anchor w -side top -padx 3 -pady 1
}

pack forget $w.contrastest $w.sscleanup
if { $fmri(level) == 1 } {
    pack $w.contrastest -in $fmri(miscf) -anchor w -side top -padx 3 -pady 1
} else {
    pack $w.sscleanup -in $fmri(miscf) -anchor w -side top -padx 5 -pady 1
}    

if { $fmri(level)==2 && $fmri(analysis)==4 } {
    set fmri(inputtype) 1
}

pack forget $w.inputtype $w.datamain $w.outputdir
if { $fmri(analysis) != 0 && $fmri(analysis) != 4 } {
    pack $w.outputdir -in $fmri(dataf) -anchor w -side top -pady 10
    if { $fmri(level) == 1  } {
	pack $w.datamain  -in $fmri(dataf) -anchor w -side top -pady 3
    } else {
	pack $w.inputtype -in $fmri(dataf) -before $w.multiple -anchor w -side top -pady 3
    }
}

#}}}

    feat5:updatestats $w 0 -1
    feat5:updatepoststats $w 0
    feat5:updatereg $w

    $w.nb raise data
}

#}}}
#{{{ feat5:updateselect

proc feat5:updateselect { w dummy } {
    global fmri feat_files

    $w.multiple.number configure -min 1
    if { $fmri(level)==2 && $fmri(analysis)!=4 } {
	$w.multiple.number configure -min 2
	if { $fmri(multiple) < 2 } {
	    set fmri(multiple) 2
	}
    }

    if { $fmri(level) == 1 && $fmri(analysis) != 4 && $fmri(analysis) != 0 } {
	$w.multiple.setup configure -text "Select 4D data"
    } else {
	if { $fmri(level) == 1 || $fmri(inputtype) == 1 } {
	    if { $fmri(multiple) == 1 } {
		$w.multiple.setup configure -text "Select FEAT directory"
	    } else {
		$w.multiple.setup configure -text "Select FEAT directories"
	    }
	} else {
		$w.multiple.setup configure -text "Select cope images"
	}
    }

    pack forget $fmri(unwarpff).unwarpmultiple $fmri(unwarpff).unwarpsingle $fmri(unwarpff).unwarpmagmultiple $fmri(unwarpff).unwarpmagsingle \
	    $fmri(initial_highresf).initial_highresmultiple $fmri(initial_highresf).initial_highressingle \
	    $fmri(highresf).highresmultiple $fmri(highresf).highressingle

    if { ! $fmri(relative_yn) } {
	if { $fmri(multiple) < 2 } {
	    pack $fmri(unwarpff).unwarpsingle $fmri(unwarpff).unwarpmagsingle -in $fmri(unwarpff) -before $fmri(unwarpff).opts1 -anchor w -side top -pady 2 -padx 3
	    pack $fmri(initial_highresf).initial_highressingle -in $fmri(initial_highresf) -before $fmri(initial_highresf).opts -anchor w -side top -pady 2 -padx 3
	    pack $fmri(highresf).highressingle -in $fmri(highresf) -before $fmri(highresf).opts -anchor w -side top -pady 2 -padx 3
	} else {
	    pack $fmri(unwarpff).unwarpmultiple $fmri(unwarpff).unwarpmagmultiple -in $fmri(unwarpff) -before $fmri(unwarpff).opts1 -anchor w -side top -pady 2 -padx 3
	    pack $fmri(initial_highresf).initial_highresmultiple -in $fmri(initial_highresf) -before $fmri(initial_highresf).opts -anchor w -side top -pady 2 -padx 3
	    pack $fmri(highresf).highresmultiple -in $fmri(highresf) -before $fmri(highresf).opts -anchor w -side top -pady 2 -padx 3
	}
    }

    if { [ winfo exists $w.copeinputs ] } {
	destroy $w.copeinputs
    }
    if { $fmri(level) == 2 &&
	 $fmri(inputtype) == 1 &&
	 $fmri(analysis) != 4 &&
	 [ info exists feat_files(1) ] &&
	 [ file exists $feat_files(1) ] } {
	#{{{ input cope select

frame $w.copeinputs
pack $w.copeinputs -in $fmri(dataf) -before $w.outputdir -side top -anchor w -padx 4 -pady 4

set w0 $w.copeinputs
frame $w0.f
pack $w0.f -expand yes -fill both -in $w0 -side top
canvas $w0.f.viewport -xscrollcommand "$w0.xsbar set" -borderwidth 0
scrollbar $w0.xsbar -command "$w0.f.viewport xview" -orient horizontal
frame $w0.f.viewport.f
$w0.f.viewport create window 0 0 -anchor nw -window $w0.f.viewport.f
bind $w0.f.viewport.f <Configure> "feat5:scrollform_resize $w0 $w0.f.viewport"
pack $w0.f.viewport -side left -fill both -expand true -in $w0.f
set v $w0.f.viewport.f

set statslist [ lsort -dictionary [ imglob -oneperimage $feat_files(1)/stats/cope*.* ] ]
set fmri(ncopeinputs) [ llength $statslist ]
    
if { $fmri(ncopeinputs) < 1 } {
    MxPause "Warning: the first selected FEAT directory contains no stats/cope images"
    return 1
}

label $v.0 -text "Use lower-level copes: "
grid $v.0 -in $v -column 0 -row 0

for { set nci 1 } { $nci <=  $fmri(ncopeinputs) } { incr nci 1 } {   
    if { ! [ info exists fmri(copeinput.$nci) ] } {
	set fmri(copeinput.$nci) 1
    }
    checkbutton $v.$nci -variable fmri(copeinput.$nci) -text "$nci "
    grid $v.$nci -in $v -column $nci -row 0
}

$w.bhelp bind $w.copeinputs -msg "The higher-level FEAT analysis will be run separately for each
lower-level contrast. You can tell FEAT to ignore certain lower-level
contrasts by turning off the appropriate button."

#}}}
    }

}

#}}}
#{{{ feat5:updatehelp

proc feat5:updatehelp { w } {
    global fmri

    if { $fmri(help_yn) == 1 } {
	$w.bhelp configure -state balloon
    } else {
	$w.bhelp configure -state none
    }
}

#}}}
#{{{ feat5:updateperfusion

proc feat5:updateperfusion { w } {
    global fmri

    set fmri(templp_yn) 0

    if { ! $fmri(perfsub_yn) } {
        set fmri(prewhiten_yn) 1
        pack forget $w.temp.tcmenu
    } else {
        set fmri(prewhiten_yn) 0
        pack $w.temp.tcmenu -in $w.temp -after $w.temp.ps_yn -side top -side left -padx 5
    }

    MxPause "Warning - you have changed the \"Perfusion subtraction\" setting, which may have changed the prewhitening option."

}

#}}}
#{{{ feat5:updateprestats

proc feat5:updateprestats { w dummy } {
    global fmri

    if { $fmri(st) < 3 || $fmri(st) > 4 } {
	pack forget $w.st_file
    } else {
	pack $w.st_file -in $fmri(stf) -side left -padx 5
    }

    if { $fmri(regunwarp_yn) } {
	pack forget $fmri(unwarpf).label 
	pack $fmri(unwarpf).lf -in $fmri(unwarpf) -side left
    } else {
	pack forget $fmri(unwarpf).lf 
	pack $fmri(unwarpf).label -in $fmri(unwarpf) -side left -before $fmri(unwarpf).yn
    }

}

#}}}
#{{{ feat5:updatestats

proc feat5:updatestats { w process dummy } {
    global fmri

    if { $process } {
	if { [ info exists fmri(w_model) ] && [ winfo exists $fmri(w_model) ] } {
	    destroy $fmri(w_model)
	}
	feat5:setup_model_vars_simple $w
	if { [ feat5:setup_model_preview $w ] } {
	    return 1
	}
    }

    if { $fmri(infeat) } {
	pack forget $w.prewhiten $w.motionevs $w.wizard $w.model $w.mixed
	if { $fmri(level) == 1 } {
	    pack $w.prewhiten $w.motionevs $w.wizard $w.model -in $fmri(statsf) -anchor w -side top -padx 5 -pady 3
	} else {
	    pack $w.mixed     $w.wizard $w.model -in $fmri(statsf) -anchor w -side top -padx 5 -pady 3
	}
    } else {
	glm:updatelevel $w 0
    }
}

#}}}
#{{{ feat5:updatepoststats

proc feat5:updatepoststats { w dummy } {
    global fmri

    pack forget $w.modelcon $w.uncorrected $w.voxel $w.z_thresh $w.prob_thresh $w.conmask $w.render $w.bgimage $w.zmin $w.zmax

    if { ! $fmri(stats_yn) } {
	pack $w.modelcon -in $fmri(poststatsf) -before $fmri(poststatsf).threshmask -side top -anchor w -padx 5 -pady 5
    }

    if { $fmri(thresh) } {
	if { $fmri(thresh) == 1 } {
	    pack $w.uncorrected -in $fmri(lfthresh) -side left
	} elseif { $fmri(thresh) == 2 } {
	    pack $w.voxel -in $fmri(lfthresh) -side left
	} else {
	    pack $w.z_thresh $w.prob_thresh -in $fmri(lfthresh) -side left -anchor w -pady 2
	}

	pack $w.conmask -in $fmri(poststatsf) -after $w.thresh -side top -anchor w -pady 5 -padx 5
	pack $w.render  -in $fmri(poststatsf) -after $w.conmask -side top -anchor w -pady 5

	pack $fmri(lfrenderingtop) -in $fmri(lfrendering) -anchor w -padx 3 -pady 3

	if { $fmri(level) > 1 } {
	    pack $w.bgimage -in $fmri(lfrendering) -before $fmri(lfrenderingtop) -anchor w -padx 3 -pady 3
	}
	
	if { $fmri(zdisplay) } {
	    pack $w.zmin $w.zmax -in $fmri(lfrenderingtop) -after $w.zmaxmenu -side left -anchor w -padx 3 -pady 5
	}
    }

}

#}}}
#{{{ feat5:updatereg

proc feat5:updatereg_hr_init { w } {
    global fmri

    if { $fmri(reginitial_highres_yn) == 1 } {
	set fmri(reghighres_yn) 1
    }

    feat5:updatereg $w
}

proc feat5:updatereg_hr { w } {
    global fmri

    if { $fmri(reghighres_yn) == 0 } {
	set fmri(reginitial_highres_yn) 0
    }

    feat5:updatereg $w
}

proc feat5:updatereg { w } {
    global fmri

    if { $fmri(regreduce_dof) > 0 } {
	set fmri(reginitial_highres_dof) 12
	set fmri(reghighres_dof) 12
	set fmri(regstandard_dof) 12

	set thedof 7
	if { $fmri(regreduce_dof) == 2 } {
	    set thedof 3
	}

	if { $fmri(reginitial_highres_yn) } {
	    set fmri(reginitial_highres_dof) $thedof
	} elseif { $fmri(reghighres_yn) } {
	    set fmri(reghighres_dof) $thedof
	} elseif { $fmri(regstandard_yn) } {
	    set fmri(regstandard_dof) $thedof
	}
    }

    if { $fmri(reginitial_highres_yn) } {
	pack forget $fmri(regf).initial_highres.label 
	pack $fmri(regf).initial_highres.lf -in $fmri(regf).initial_highres -side left
    } else {
	pack forget $fmri(regf).initial_highres.lf 
	pack $fmri(regf).initial_highres.label -in $fmri(regf).initial_highres -side left
    }

    if { $fmri(reghighres_yn) } {
	pack forget $fmri(regf).highres.label 
	pack $fmri(regf).highres.lf -in $fmri(regf).highres -side left
    } else {
	pack forget $fmri(regf).highres.lf 
	pack $fmri(regf).highres.label -in $fmri(regf).highres -side left
    }

    if { $fmri(regstandard_yn) } {
	pack forget $fmri(regf).standard.label 
	pack $fmri(regf).standard.lf -in $fmri(regf).standard -side left
    } else {
	pack forget $fmri(regf).standard.lf 
	pack $fmri(regf).standard.label -in $fmri(regf).standard -side left
    }

    feat5:updateselect $w 0
}

#}}}
#{{{ feat5:updateimageinfo

proc feat5:updateimageinfo { w i full } {
    global FSLDIR feat_files fmri

    set thefile [ remove_ext $feat_files($i) ]
    set changed_stuff 0

    if { [ imtest $thefile ] } {
	set npts [ exec sh -c "${FSLDIR}/bin/avwnvols $thefile 2> /dev/null" ]

	if { $npts > 0 } {
	    set fmri(npts) $npts
	}

	if { $full } {

	    # set BET and FLIRT DOF according to FOV
	    set xfov [ expr abs([ exec sh -c "$FSLDIR/bin/avwval $thefile pixdim1" ] * [ exec sh -c "$FSLDIR/bin/avwval $thefile dim1" ]) ]
	    set yfov [ expr abs([ exec sh -c "$FSLDIR/bin/avwval $thefile pixdim2" ] * [ exec sh -c "$FSLDIR/bin/avwval $thefile dim2" ]) ]
	    set zfov [ expr abs([ exec sh -c "$FSLDIR/bin/avwval $thefile pixdim3" ] * [ exec sh -c "$FSLDIR/bin/avwval $thefile dim3" ]) ]
	    
	    set BETMINFOV 30
	    if { $xfov < $BETMINFOV || $yfov < $BETMINFOV || $zfov < $BETMINFOV } {
		set fmri(bet_yn) 0
		set changed_stuff 1
	    }

	    set fmri(regreduce_dof) 0

	    set REGMINFOV 120
	    if { $xfov < $REGMINFOV || $yfov < $REGMINFOV || $zfov < $REGMINFOV } {
		set fmri(regreduce_dof) 1
		feat5:updatereg $w
		set changed_stuff 1
	    }

	    set REGTRANSMINFOV 60
	    if { $xfov < $REGTRANSMINFOV || $yfov < $REGTRANSMINFOV || $zfov < $REGTRANSMINFOV } {
		set fmri(regreduce_dof) 2
		feat5:updatereg $w
		set changed_stuff 1
	    }

	    if { $changed_stuff } {
		MxPause "Warning - have auto-set BET preprocessing option and/or registration DoF on the basis of image fields-of-view; check settings."
	    }
	}
    }

    if { $fmri(level) > 1 && $fmri(analysis) != 4 } {
	set fmri(npts) $fmri(multiple)
    }
}

#}}}
#{{{ feat5:write

proc feat5:write { w feat_model write_image_filenames exitoncheckfail filename } {
    global fmri FSLDIR feat_files unwarp_files unwarp_files_mag initial_highres_files highres_files

    if { $fmri(level) == 1 && $fmri(analysis) > 0 && $fmri(con_mode) == "orig" && [ feat5:setup_model_update_contrasts_mode $w 0 ] == -1 } {
	return -1
    }

    set filename [ file rootname $filename ]

    if { $fmri(level) > 1 } {
	set fmri(npts) $fmri(multiple)
	set fmri(ndelete) 0
    }

    set channel [ open ${filename}.fsf "w" ]

    #{{{ basic variables

    puts $channel "
# FEAT version number
set fmri(version) $fmri(version)

# Analysis level
# 1 : First-level analysis
# 2 : Higher-level analysis
set fmri(level) $fmri(level)

# Which stages to run
# 0 : No first-level analysis (registration and/or group stats only)
# 7 : Full first-level analysis
# 1 : Pre-Stats
# 3 : Pre-Stats + Stats
# 2 :             Stats
# 6 :             Stats + Post-stats
# 4 :                     Post-stats
set fmri(analysis) $fmri(analysis)

# Delay before starting (hours)
set fmri(delay) $fmri(delay)

# Use relative filenames
set fmri(relative_yn) $fmri(relative_yn)

# Balloon help
set fmri(help_yn) $fmri(help_yn)

# Run Featwatcher
set fmri(featwatcher_yn) $fmri(featwatcher_yn)

# Cleanup first-level standard-space images
set fmri(sscleanup_yn) $fmri(sscleanup_yn)

# Output directory
set fmri(outputdir) \"$fmri(outputdir)\"

# TR(s)
set fmri(tr) $fmri(tr)

# Total volumes
set fmri(npts) $fmri(npts)

# Delete volumes
set fmri(ndelete) $fmri(ndelete)

# Perfusion tag/control order
set fmri(tagfirst) $fmri(tagfirst)

# Number of first-level analyses
set fmri(multiple) $fmri(multiple)

# Higher-level input type
# 1 : Inputs are lower-level FEAT directories
# 2 : Inputs are cope images from FEAT directories
set fmri(inputtype) $fmri(inputtype)

# Carry out pre-stats processing?
set fmri(filtering_yn) $fmri(filtering_yn)

# Brain/background threshold, %
set fmri(brain_thresh) $fmri(brain_thresh)

# Critical z for design efficiency calculation
set fmri(critical_z) $fmri(critical_z)

# Noise level
set fmri(noise) $fmri(noise)

# Noise AR(1)
set fmri(noisear) $fmri(noisear)

# Post-stats-only directory copying
# 0 : Overwrite original post-stats results
# 1 : Copy original FEAT directory for new Contrasts, Thresholding, Rendering
set fmri(newdir_yn) $fmri(newdir_yn)

# Motion correction
# 0 : None
# 1 : MCFLIRT
set fmri(mc) $fmri(mc)

# Spin-history (currently obsolete)
set fmri(sh_yn) $fmri(sh_yn)

# B0 fieldmap unwarping?
set fmri(regunwarp_yn) $fmri(regunwarp_yn)

# EPI dwell time (ms)
set fmri(dwell) $fmri(dwell)

# EPI TE (ms)
set fmri(te) $fmri(te)

# % Signal loss threshold
set fmri(signallossthresh) $fmri(signallossthresh)

# Unwarp direction
set fmri(unwarp_dir) $fmri(unwarp_dir)

# Slice timing correction
# 0 : None
# 1 : Regular up (0, 1, 2, 3, ...)
# 2 : Regular down
# 3 : Use slice order file
# 4 : Use slice timings file
# 5 : Interleaved (0, 2, 4 ... 1, 3, 5 ... )
set fmri(st) $fmri(st)

# Slice timings file
set fmri(st_file) \"$fmri(st_file)\"

# BET brain extraction
set fmri(bet_yn) $fmri(bet_yn)

# Spatial smoothing FWHM (mm)
set fmri(smooth) $fmri(smooth)

# Intensity normalization
set fmri(norm_yn) $fmri(norm_yn)

# Perfusion subtraction
set fmri(perfsub_yn) $fmri(perfsub_yn)

# Highpass temporal filtering
set fmri(temphp_yn) $fmri(temphp_yn)

# Lowpass temporal filtering
set fmri(templp_yn) $fmri(templp_yn)

# MELODIC ICA data exploration
set fmri(melodic_yn) $fmri(melodic_yn)

# Carry out main stats?
set fmri(stats_yn) $fmri(stats_yn)

# Carry out prewhitening?
set fmri(prewhiten_yn) $fmri(prewhiten_yn)

# Add motion parameters to model
# 0 : No
# 1 : Yes, orthogonalise rest of model wrt motion
# 2 : Yes, orthogonalise motion wrt rest of model
# 3 : Yes, don't orthogonalise
set fmri(motionevs) $fmri(motionevs)

# Higher-level modelling
# 3 : Fixed effects
# 0 : Mixed Effects: Simple OLS
# 2 : Mixed Effects: FLAME 1
# 1 : Mixed Effects: FLAME 1+2
set fmri(mixed_yn) $fmri(mixed_yn)

# Number of EVs
set fmri(evs_orig) $fmri(evs_orig)
set fmri(evs_real) $fmri(evs_real)

# Number of contrasts
set fmri(ncon_orig) $fmri(ncon_orig)
set fmri(ncon_real) $fmri(ncon_real)

# Number of F-tests
set fmri(nftests_orig) $fmri(nftests_orig)
set fmri(nftests_real) $fmri(nftests_real)

# Add constant column to design matrix? (obsolete)
set fmri(constcol) $fmri(constcol)

# Carry out post-stats steps?
set fmri(poststats_yn) $fmri(poststats_yn)

# Pre-threshold masking?
set fmri(threshmask) \"$fmri(threshmask)\"

# Thresholding
# 0 : None
# 1 : Uncorrected
# 2 : Voxel
# 3 : Cluster
set fmri(thresh) $fmri(thresh)

# P threshold
set fmri(prob_thresh) $fmri(prob_thresh)

# Z threshold
set fmri(z_thresh) $fmri(z_thresh)

# Z min/max for colour rendering
# 0 : Use actual Z min/max
# 1 : Use preset Z min/max
set fmri(zdisplay) $fmri(zdisplay)

# Z min in colour rendering
set fmri(zmin) $fmri(zmin)

# Z max in colour rendering
set fmri(zmax) $fmri(zmax)

# Colour rendering type
# 0 : Solid blobs
# 1 : Transparent blobs
set fmri(rendertype) $fmri(rendertype)

# Background image for higher-level stats overlays
# 1 : Mean highres
# 2 : First highres
# 3 : Mean functional
# 4 : First functional
# 5 : Standard space template
set fmri(bgimage) $fmri(bgimage)

# Registration?
set fmri(reg_yn) $fmri(reg_yn)

# Registration to initial structural
set fmri(reginitial_highres_yn) $fmri(reginitial_highres_yn)

# Search space for registration to initial structural
# 0   : No search
# 90  : Normal search
# 180 : Full search
set fmri(reginitial_highres_search) $fmri(reginitial_highres_search)

# Degrees of Freedom for registration to initial structural
set fmri(reginitial_highres_dof) $fmri(reginitial_highres_dof)

# Do nonlinear registration to initial structural?
set fmri(reginitial_highres_nonlinear_yn) $fmri(reginitial_highres_nonlinear_yn)

# Registration to main structural
set fmri(reghighres_yn) $fmri(reghighres_yn)

# Search space for registration to main structural
# 0   : No search
# 90  : Normal search
# 180 : Full search
set fmri(reghighres_search) $fmri(reghighres_search)

# Degrees of Freedom for registration to main structural
set fmri(reghighres_dof) $fmri(reghighres_dof)

# Do nonlinear registration to main structural?
set fmri(reghighres_nonlinear_yn) $fmri(reghighres_nonlinear_yn)

# Registration to standard image?
set fmri(regstandard_yn) $fmri(regstandard_yn)

# Standard image
set fmri(regstandard) \"$fmri(regstandard)\"

# Search space for registration to standard space
# 0   : No search
# 90  : Normal search
# 180 : Full search
set fmri(regstandard_search) $fmri(regstandard_search)

# Degrees of Freedom for registration to standard space
set fmri(regstandard_dof) $fmri(regstandard_dof)

# Do nonlinear registration to standard space?
set fmri(regstandard_nonlinear_yn) $fmri(regstandard_nonlinear_yn)

# High pass filter cutoff
set fmri(paradigm_hp) $fmri(paradigm_hp)"

#}}}
    #{{{ input and highres filenames

puts $channel "
# Number of lower-level copes feeding into higher-level analysis
set fmri(ncopeinputs) $fmri(ncopeinputs)"
for { set nci 1 } { $nci <=  $fmri(ncopeinputs) } { incr nci 1 } {   
	puts $channel "
# Use lower-level cope $nci for higher-level analysis
set fmri(copeinput.$nci) $fmri(copeinput.$nci)"
}

if { $write_image_filenames } {

if { $fmri(multiple) > 0 } {
    if { $w != -1 } {
	if { [ feat5:multiple_check $w 0 0 0 0 ] && $exitoncheckfail } {
	    return 1
	}
    }
    for { set i 1 } { $i <= $fmri(multiple) } { incr i 1 } {
	puts $channel "
# 4D AVW data or FEAT directory ($i)
set feat_files($i) \"$feat_files($i)\""
    }
}


if { $fmri(reg_yn) && ( $fmri(regunwarp_yn) || $fmri(reginitial_highres_yn) || $fmri(reghighres_yn) ) } {

    if { $fmri(multiple) < 2 } {
	set nhighres 1
    } else {
	set nhighres $fmri(multiple)
    }

    if { $fmri(relative_yn) } {
	for { set i 1 } { $i <= $nhighres } { incr i 1 } {
	    set unwarp_files($i) [ file dirname $feat_files($i) ]/blahstruct_brain.hdr
	    set initial_highres_files($i) [ file dirname $feat_files($i) ]/struct_brain.hdr
	    set highres_files($i) [ file dirname $feat_files($i) ]/blah2struct_brain.hdr
	}
    }

    if { $w != -1 } {
	if { $fmri(regunwarp_yn) } {
	    if { [ feat5:multiple_check $w 1 0 0 0 ] && $exitoncheckfail } {
		return 1
	    }
	    for { set i 1 } { $i <= $nhighres } { incr i 1 } {
		puts $channel "
# B0 unwarp input image for analysis $i
set unwarp_files($i) \"$unwarp_files($i)\""
            }
	    if { [ feat5:multiple_check $w 2 0 0 0 ] && $exitoncheckfail } {
		return 1
	    }
	    for { set i 1 } { $i <= $nhighres } { incr i 1 } {
		puts $channel "
# B0 unwarp mag input image for analysis $i
set unwarp_files_mag($i) \"$unwarp_files_mag($i)\""
            }
	}
	if { $fmri(reginitial_highres_yn) } {
	    if { [ feat5:multiple_check $w 3 0 0 0 ] && $exitoncheckfail } {
		return 1
	    }
	    for { set i 1 } { $i <= $nhighres } { incr i 1 } {
		puts $channel "
# Session's structural image for analysis $i
set initial_highres_files($i) \"$initial_highres_files($i)\""
            }
	}
	if { $fmri(reghighres_yn) } {
	    if { [ feat5:multiple_check $w 4 0 0 0 ] && $exitoncheckfail } {
		return 1
	    }
	    for { set i 1 } { $i <= $nhighres } { incr i 1 } {
		puts $channel "
# Subject's structural image for analysis $i
set highres_files($i) \"$highres_files($i)\""
            }
	}
    }
}

} else { 
    if { $fmri(npts) < 1 } {
	MxPause "Please either select input data or set the number of total volumes before attempting to view the design."
	if { $exitoncheckfail } {
	    return -1
	}
    }
}

#}}}
    #{{{ EVs

	for { set i 1 } { $i <= $fmri(evs_orig) } { incr i 1 } {

	    puts $channel "
# EV $i title
set fmri(evtitle$i) \"$fmri(evtitle$i)\""

	    if { [ info exists fmri(shape$i) ] } {

		puts $channel "
# Basic waveform shape (EV $i)
# 0 : Square
# 1 : Sinusoid
# 2 : Custom (1 entry per volume)
# 3 : Custom (3 column format)
# 4 : Interaction
# 10 : Empty (all zeros)
set fmri(shape$i) $fmri(shape$i)

# Convolution (EV $i)
# 0 : None
# 1 : Gaussian
# 2 : Gamma
# 3 : Double-Gamma HRF
# 4 : Gamma basis functions
# 5 : Sine basis functions
# 6 : FIR basis functions
set fmri(convolve$i) $fmri(convolve$i)

# Convolve phase (EV $i)
set fmri(convolve_phase$i) $fmri(convolve_phase$i)

# Apply temporal filtering (EV $i)
set fmri(tempfilt_yn$i) $fmri(tempfilt_yn$i)

# Add temporal derivative (EV $i)
set fmri(deriv_yn$i) $fmri(deriv_yn$i)"

		switch $fmri(shape$i) {
		    0 { 
			puts $channel "
# Skip (EV $i)
set fmri(skip$i) $fmri(skip$i)

# Off (EV $i)
set fmri(off$i) $fmri(off$i)

# On (EV $i)
set fmri(on$i) $fmri(on$i)

# Phase (EV $i)
set fmri(phase$i) $fmri(phase$i)

# Stop (EV $i)
set fmri(stop$i) $fmri(stop$i)"
		    }
		    1 { 
			puts $channel "
# Skip (EV $i)
set fmri(skip$i) $fmri(skip$i)

# Period (EV $i)
set fmri(period$i) $fmri(period$i)

# Phase (EV $i)
set fmri(phase$i) $fmri(phase$i)

# Sinusoid harmonics (EV $i)
set fmri(nharmonics$i) $fmri(nharmonics$i) 

# Stop (EV $i)
set fmri(stop$i) $fmri(stop$i)"
		    }
		    2 { puts $channel "
# Custom EV file (EV $i)
set fmri(custom$i) \"$fmri(custom$i)\"" }
		    3 { puts $channel "
# Custom EV file (EV $i)
set fmri(custom$i) \"$fmri(custom$i)\"" }
		    4 {
			for { set j 1 } { $j < $i } { incr j 1 } {
			    puts $channel "
# Interactions (EV $i with EV $j)
set fmri(interactions${i}.$j) $fmri(interactions${i}.$j)

# Demean before using in interactions (EV $i with EV $j)
set fmri(interactionsd${i}.$j) $fmri(interactionsd${i}.$j)"
			}
		    }
		}
		
		switch $fmri(convolve$i) {
		    0 { }
		    1 { 
			puts $channel "
# Gauss sigma (EV $i)
set fmri(gausssigma$i) $fmri(gausssigma$i)

# Gauss delay (EV $i)
set fmri(gaussdelay$i) $fmri(gaussdelay$i)"
		    }
		    2 {
			puts $channel "
# Gamma sigma (EV $i)
set fmri(gammasigma$i) $fmri(gammasigma$i)

# Gamma delay (EV $i)
set fmri(gammadelay$i) $fmri(gammadelay$i)"
		    }
		    3 { }
		    4 {
			puts $channel "
# Gamma basis functions number (EV $i)
set fmri(basisfnum$i) $fmri(basisfnum$i)

# Gamma basis functions window(s) (EV $i)
set fmri(basisfwidth$i) $fmri(basisfwidth$i)"
                    }
		    5 {
			puts $channel "
# Sine basis functions number (EV $i)
set fmri(basisfnum$i) $fmri(basisfnum$i)

# Sine basis functions window(s) (EV $i)
set fmri(basisfwidth$i) $fmri(basisfwidth$i)"
                    }
		    6 {
			puts $channel "
# FIR basis functions number (EV $i)
set fmri(basisfnum$i) $fmri(basisfnum$i)

# FIR basis functions window(s) (EV $i)
set fmri(basisfwidth$i) $fmri(basisfwidth$i)"
                    }
		    7 {
			puts $channel "
# FIR basis functions number (EV $i)
set fmri(basisfnum$i) $fmri(basisfnum$i)

# Optimal/custom HRF convolution file (EV $i)
set fmri(bfcustom$i) \"$fmri(bfcustom$i)\""
                    }
		}

		for { set j 0 } { $j <= $fmri(evs_orig) } { incr j 1 } {
		    puts $channel "
# Orthogonalise EV $i wrt EV $j
set fmri(ortho${i}.$j) $fmri(ortho${i}.$j)"
		}

		if { $fmri(level) > 1 } {
		    for { set j 1 } { $j <= $fmri(npts) } { incr j 1 } {
			puts $channel "
# Higher-level EV value for EV $i and input $j
set fmri(evg${j}.$i) $fmri(evg${j}.$i)"
		    }
		}
		
	    }
	}

	if { $fmri(level) > 1 } {
	    if { [ info exists fmri(level2orth) ] } {
		puts $channel "
# Setup orthogonalisation at higher level? 
set fmri(level2orth) $fmri(level2orth)"
	    }

            for { set j 1 } { $j <= $fmri(multiple) } { incr j 1 } {
		puts $channel "
# Group membership for input $j
set fmri(groupmem.$j) $fmri(groupmem.$j)"
            }
	}

#}}}
    #{{{ contrasts & F-tests

puts $channel "
# Contrast & F-tests mode
# real : control real EVs
# orig : control original EVs
set fmri(con_mode_old) $fmri(con_mode)
set fmri(con_mode) $fmri(con_mode)"

if { $fmri(level) == 1 } {
    set modes "real orig"
} else {
    set modes real
}

foreach cm $modes {

    for { set i 1 } { $i <= $fmri(ncon_${cm}) } { incr i 1 } {

	puts $channel "
# Display images for contrast_${cm} $i
set fmri(conpic_${cm}.$i) $fmri(conpic_${cm}.$i)

# Title for contrast_${cm} $i
set fmri(conname_${cm}.$i) \"$fmri(conname_${cm}.$i)\""

        for { set j 1 } { $j <= $fmri(evs_${cm}) } { incr j 1 } {
	    puts $channel "
# Real contrast_${cm} vector $i element $j
set fmri(con_${cm}${i}.${j}) $fmri(con_${cm}${i}.${j})"
        }

	for { set j 1 } { $j <= $fmri(nftests_${cm}) } { incr j 1 } {
	    puts $channel "
# F-test $j element $i
set fmri(ftest_${cm}${j}.${i}) $fmri(ftest_${cm}${j}.${i})"
        }

    }

}

puts $channel "
# Contrast masking - use >0 instead of thresholding?
set fmri(conmask_zerothresh_yn) $fmri(conmask_zerothresh_yn)"

set total [ expr $fmri(ncon_real) + $fmri(nftests_real) ]
set fmri(conmask1_1) 0
for { set c 1 } { $c <= $total } { incr c } {
    for { set C 1 } { $C <= $total } { incr C } {
	if { $C != $c } {
	    if { ! [ info exists fmri(conmask${c}_${C}) ] } {
		set fmri(conmask${c}_${C}) 0
	    }
	    puts $channel "
# Mask real contrast/F-test $c with real contrast/F-test $C?
set fmri(conmask${c}_${C}) $fmri(conmask${c}_${C})"

	    if { $fmri(conmask${c}_${C}) } {
		set fmri(conmask1_1) 1
	    }
	}
    }
}

puts $channel "
# Do contrast masking at all?
set fmri(conmask1_1) $fmri(conmask1_1)"

#}}}
    #{{{ non-GUI options

    puts $channel "
##########################################################
# Now options that don't appear in the GUI

# Alternative example_func image (not derived from input 4D dataset)
set fmri(alternative_example_func) \"$fmri(alternative_example_func)\"

# Delete GLM residuals?
set fmri(cleanup_residuals_yn) $fmri(cleanup_residuals_yn)

# Initial structural space registration initialisation transform
set fmri(init_initial_highres) \"$fmri(init_initial_highres)\"

# Structural space registration initialisation transform
set fmri(init_highres) \"$fmri(init_highres)\"

# Standard space registration initialisation transform
set fmri(init_standard) \"$fmri(init_standard)\""

#}}}

    close $channel
    
    if { $w != -1 } {
	set result 0
	if { $feat_model } {

	    set result [ catch { exec sh -c "${FSLDIR}/bin/feat_model $filename" } ErrMsg ]
	    if {$result != 0 || [ string length $ErrMsg ] > 0 } {
		    MxPause "Problem with processing the model: $ErrMsg"
	    }
	}

	return $result
    }
}

#}}}
#{{{ feat5:load

proc feat5:load { w full filename } {
    global fmri feat_files unwarp_files unwarp_files_mag initial_highres_files highres_files

    set FEATVERSION $fmri(version)

    set version [ exec sh -c "grep 'fmri(version)' $filename | awk '{ print \$3 }'" ]

    if { $version > 4.99 } {

	if { $w != -1 } {
	    set level    $fmri(level)
	    set analysis $fmri(analysis)
	    set multiple $fmri(multiple)
	    set reg_yn   $fmri(reg_yn)
	    for { set i 1 } { $i <= $multiple } { incr i 1 } {
		if { [ info exists feat_files($i) ] } {
		    set lfeat_files($i) $feat_files($i)
		}
	    }
	}

	source ${filename}

	if { $w != -1 && $fmri(analysis) != 4 && $fmri(analysis) != 0 } {
	    feat5:updateimageinfo $w 1 0
	}

	if { ! $full } {
	    if { $level != $fmri(level) } {
		feat5:updatelevel $w 0
		set fmri(level) $level
	    }

	    set fmri(npts) [ expr $fmri(npts) - $fmri(ndelete) ]
	    set fmri(ndelete) 0
	    set fmri(reg_yn)   $reg_yn
	    set fmri(analysis) $analysis
	    set fmri(multiple) $multiple
	    for { set i 1 } { $i <= $multiple } { incr i 1 } {
		set feat_files($i) $lfeat_files($i)
	    }
	}

	if { $w != -1 } {
	    feat5:updateanalysis $w 0
	    feat5:updatehelp $w
	    feat5:updateprestats $w 0
	}

    } else {
	MxPause "FEAT setup file is too old to load - sorry!"
    }
    
    set fmri(version) $FEATVERSION
}

#}}}
#{{{ feat5:wizard

proc feat5:wizard { w } {
    global FSLDIR fmri

    #{{{ setup window

    set count 0
    set w0 ".dialog[incr count]"
    while { [ winfo exists $w0 ] } {
        set w0 ".dialog[incr count]"
    }

    toplevel $w0

    wm title $w0 "Model setup wizard"
    wm iconname $w0 "wizard"
    wm iconbitmap $w0 @${FSLDIR}/tcl/fmrib.xbm

    frame $w0.f
    pack $w0.f

#}}}

    if { $fmri(level) == 1 } {
	#{{{ choose paradigm type

frame $w0.f.paradigm_type

set fmri(wizard_type) 1	

radiobutton $w0.f.paradigm_type.jab -text "rArA..." -value 1 -variable fmri(wizard_type) -command "feat5:update_wizard $w0"
radiobutton $w0.f.paradigm_type.jabac -text "rArBrArB..." -value 2 -variable fmri(wizard_type) -command "feat5:update_wizard $w0"
radiobutton $w0.f.paradigm_type.jperfab -text "perfusion rArA..." -value 3 -variable fmri(wizard_type) -command "feat5:update_wizard $w0"

$w.bhelp bind $w0.f.paradigm_type -msg "Choose whether to setup rArA... or rArBrArB... designs (regular block
or single-event). The r blocks will normally be rest (or control)
conditions.

The \"perfusion rArA...\" option sets up the full model for a simple
perfusion experiment, setting up a constant-height control-tag EV, an
average BOLD component EV and the interaction EV, which represents the
control-tag functional modulation."

pack $w0.f.paradigm_type.jab $w0.f.paradigm_type.jabac $w0.f.paradigm_type.jperfab -in $w0.f.paradigm_type -side top -padx 5 -side left

#}}}
	#{{{ frame ab

frame $w0.f.ab

set fmri(r_count) 30
tixControl $w0.f.ab.r_count -label "r (rest) period (s)" \
	-variable fmri(r_count) -step 1 -min 0 -selectmode immediate -options { entry.width 4 }

set fmri(a_count) 30
tixControl $w0.f.ab.a_count -label "A period (s)" \
	-variable fmri(a_count) -step 1 -min 0 -selectmode immediate -options { entry.width 4 }

set fmri(b_count) 30
tixControl $w0.f.ab.b_count -label "B period (s)" \
	-variable fmri(b_count) -step 1 -min 0 -selectmode immediate -options { entry.width 4 }

pack $w0.f.ab.r_count $w0.f.ab.a_count -in $w0.f.ab -side left -padx 5 -pady 3 -anchor n

$w.bhelp bind $w0.f.ab.r_count -msg "The r period (seconds) is normally the rest period, i.e., no
stimulation was applied during this period."

$w.bhelp bind $w0.f.ab.a_count -msg "The A period (seconds) is normally the activation period, i.e.,
stimulation was applied during this period."

$w.bhelp bind $w0.f.ab.b_count -msg "The B period (seconds) is normally the period associated with a second
type of activation, i.e., stimulation type 2 was applied during this
period."

#}}}
    } else {
	#{{{ choose paradigm type

frame $w0.f.paradigm_type

set fmri(wizard_type) 1	

radiobutton $w0.f.paradigm_type.onegroup -text "single group average" -value 1 -variable fmri(wizard_type) -command "feat5:update_wizard $w0"
radiobutton $w0.f.paradigm_type.twogroupunpaired -text "two groups, unpaired" -value 2 -variable fmri(wizard_type) -command "feat5:update_wizard $w0"
radiobutton $w0.f.paradigm_type.twogrouppaired -text "two groups, paired" -value 3 -variable fmri(wizard_type) -command "feat5:update_wizard $w0"

$w.bhelp bind $w0.f.paradigm_type -msg "Select from the different higher-level designs. In the case of the
unpaired two-group test, set the number of subjects in the first
group. Note that in the case of the paired two-group test, your
subjects should be entered: first all subjects for the first
condition, then all the subjects (in the same order) for the second
condition."

pack $w0.f.paradigm_type.onegroup $w0.f.paradigm_type.twogroupunpaired $w0.f.paradigm_type.twogrouppaired -in $w0.f.paradigm_type -side top -padx 3 -anchor w

#}}}
	#{{{ frame ab

frame $w0.f.ab

tixControl $w0.f.ab.a_count -label "Number of subjects in first group" \
    -variable fmri(a_count) -step 1 -min 1 -max [ expr $fmri(multiple) - 1 ] -selectmode immediate -options { entry.width 4 }

$w.bhelp bind $w0.f.ab.a_count -msg "This is the number of subjects in the first group of subjects. The
number of subjects in the second group is calculated automatically."

#}}}
    }
    #{{{ setup button

button $w0.f.cancel -command "feat5:updatestats $w 1 -1 ; destroy $w0" -text "Process"

#}}}
 
    pack $w0.f.paradigm_type $w0.f.ab $w0.f.cancel -in $w0.f -pady 5
}

#{{{ feat5:update_wizard

proc feat5:update_wizard { w } {
    global fmri

    if { $fmri(level) == 1 } {
	if { $fmri(wizard_type) == 2 } {
	    set fmri(r_count) 30
	    set fmri(a_count) 30
	    set fmri(b_count) 30
	    pack $w.f.ab.b_count -in $w.f.ab -after $w.f.ab.a_count -side left -padx 5 -pady 3 -anchor n
	} else {
	    set fmri(r_count) 30
	    set fmri(a_count) 30
	    pack forget $w.f.ab.b_count
	}
    } else {
	pack forget $w.f.ab.a_count
	if { $fmri(wizard_type) == 2 } {
	    set fmri(a_count) 1
	    pack $w.f.ab.a_count -in $w.f.ab -side left -padx 5 -pady 3 -anchor n
	}
    }
}

#}}}

#}}}
#{{{ feat5:relative (DELETE)

proc feat5:relative { input_1 input_2 output_1 } {

    # find beginning of different substring
    for { set start 1 } { $start < [ string length $input_1 ] && $start < [ string length $output_1 ] } { incr start 1 } {
	if { [ string range $input_1  $start $start ] != [ string range $output_1 $start $start ] } {
	    break
	}
    }

    # find end of different substring
    for { set end 1 } { $end < [ string length $input_1 ] && $end < [ string length $output_1 ] } { incr end 1 } {
	if { [ string range $input_1 [ expr [ string length $input_1 ] - $end - 1 ] [ expr [ string length $input_1 ] - $end - 1 ] ] != \
		[ string range $output_1 [ expr [ string length $output_1 ] - $end - 1 ] [ expr [ string length $output_1 ] - $end - 1 ] ] } {
	    break
	}
    }

    # set substrings_1a  [ string range $input_1  0 [ expr $start - 1 ] ]
    set substrings_1bi [ string range $input_1  $start [ expr [ string length $input_1  ] - $end - 1 ] ]
    set substrings_1bo [ string range $output_1 $start [ expr [ string length $output_1 ] - $end - 1 ] ]
    # set substrings_1c  [ string range $input_1  [ expr [ string length $input_1 ] - $end ] [ expr [ string length $input_1 ] - 1 ] ]
    
    regsub $substrings_1bi $input_2 $substrings_1bo output_2

    return $output_2
}

proc feat5:relativetest { } {

feat5:relative /data/heidi/ac/ac_left.hdr /data/heidi/ac/ac_right.hdr /data/heidi/ac/ac_struct_brain.hdr
feat5:relative /data/heidi/ac/ac_left.hdr /data/heidi/at/at_left.hdr /data/heidi/ac/ac_struct_brain.hdr
feat5:relative /data/heidi/ac/ac_left.hdr /data/heidi/at/at_right.hdr /data/heidi/ac/ac_struct_brain.hdr
feat5:relative /data/heidi/ac/fmri.hdr /data/heidi/at/fmri.hdr /data/heidi/ac/struct_brain.hdr
feat5:relative /data/heidi/fmri/ac.hdr /data/heidi/fmri/at.hdr /data/heidi/struct/ac.hdr
feat5:relative /grot/s1_fmri.hdr /grot/s2_fmri.hdr /grot/s1_struct.hdr

}

#}}}
#{{{ feat5:setup_model

proc feat5:setup_model { w } {
    global FSLDIR fmri VARS PWD

    #{{{ setup window

    if { [ info exists fmri(w_model) ] } {
	if { [ winfo exists $fmri(w_model) ] } {
	    return 0
	}
    }

    set count 0
    set w0 ".wdialog[incr count]"
    while { [ winfo exists $w0 ] } {
        set w0 ".wdialog[incr count]"
    }

    set fmri(w_model) $w0

    toplevel $w0

    wm title $w0 "General Linear Model"
    wm iconname $w0 "GLM"
    wm iconbitmap $w0 @${FSLDIR}/tcl/fmrib.xbm

    frame $w0.f
    pack $w0.f -expand yes -fill both -in $w0 -side top

    canvas $w0.f.viewport -yscrollcommand "$w0.f.ysbar set" -xscrollcommand "$w0.xsbar set"
    scrollbar $w0.xsbar -command "$w0.f.viewport xview" -orient horizontal
    scrollbar $w0.f.ysbar -command "$w0.f.viewport yview" -orient vertical
    frame $w0.f.viewport.f
    $w0.f.viewport create window 0 0 -anchor nw -window $w0.f.viewport.f
    bind $w0.f.viewport.f <Configure> "feat5:scrollform_resize $w0 $w0.f.viewport"
    pack $w0.f.viewport -side left -fill both -expand true -in $w0.f

    set NB $w0.f.viewport.f.nb
    tixNoteBook $NB
    if { $fmri(analysis) != 4 } {
	$NB add evs -label "EVs"
    }
    $NB add contrasts -label "Contrasts & F-tests"
    pack $NB -in $w0.f.viewport.f

#}}}
    #{{{ setup contrasts & F-tests

set fmri(contrastsf) [ $NB subwidget contrasts ]

if { $fmri(level) == 1 } {
    tixOptionMenu $fmri(contrastsf).con_mode -label "Setup contrasts & F-tests for " -command "feat5:setup_model_update_contrasts_mode $w 1 ; feat5:setup_model_update_contrasts $w"
    $fmri(contrastsf).con_mode add command orig -label "Original EVs"
    $fmri(contrastsf).con_mode add command real -label "Real EVs"
    $fmri(contrastsf).con_mode configure -variable fmri(con_mode)
    pack $fmri(contrastsf).con_mode -in $fmri(contrastsf) -side top -anchor w -padx 5 -pady 5
    $w.bhelp bind $fmri(contrastsf).con_mode -msg "For first-level analyses, it is common for the final design matrix to
have a greater number of \"real EVs\" than the \"original\" number; for
example, when using basis functions, each \"original EV\" gives rise
to several \"real EVs\".

Therefore it is possible it many cases for you to setup contrasts and
F-tests with respect to the \"original EVs\", and FEAT will work out
for you what these will be for the final design matrix. For example, a
single \[1\] contrast on an original EV for which basis function HRF
convolutions have been chosen will result in a single \[1\] contrast
for each resulting real EV, and then an F-test across these.

In general you can switch between setting up contrasts and F-tests
with respect to \"Original EVs\" and \"Real EVs\"; though of course if
you fine-tune the contrasts for real EVs and then revert to original
EV setup some settings may be lost.

When you \"View\" the design matrix or press \"Done\" at the end of
setting up the model, an \"Original EVs\" setup will get converted to
the appropriate \"Real EVs\" settings."
}

frame $fmri(contrastsf).num
tixControl $fmri(contrastsf).num.con -label "Contrasts " \
        -variable fmri(ncon_$fmri(con_mode)) -step 1 -min 1 -command "feat5:setup_model_update_contrasts $w"
tixControl $fmri(contrastsf).num.ftests -label "   F-tests " \
        -variable fmri(nftests_$fmri(con_mode)) -step 1 -min 0 -command "feat5:setup_model_update_contrasts $w"
pack $fmri(contrastsf).num.con $fmri(contrastsf).num.ftests -in $fmri(contrastsf).num -side left -anchor n -padx 5 -pady 0
pack $fmri(contrastsf).num -in $fmri(contrastsf) -side top -anchor w -padx 0 -pady 5

feat5:setup_model_update_contrasts $w 0

$w.bhelp bind $fmri(contrastsf).num.con -msg "Each EV (explanatory variable, i.e., waveform) in the design matrix
results in a PE (parameter estimate) image.  This estimate tells you
how strongly that waveform fits the data at each voxel - the higher it
is, the better the fit. For an unblurred square wave input (which will
be scaled in the model from -0.5 to 0.5), the PE image is equivalent to
the \"mean difference image\". To convert from a PE to a t
statistic image, the PE is divided by it's standard error, which is
derived from the residual noise after the complete model has been
fit. The t image is then transformed into a Z statistic via standard
statistical transformation. As well as Z images arising from single
EVs, it is possible to combine different EVs (waveforms) - for
example, to see where one has a bigger effect than another. To do
this, one PE is subtracted from another, a combined standard error is
calculated, and a new Z image is created.

All of the above is controlled by you, by setting up contrasts. Each
output Z statistic image is generated by setting up a contrast vector;
thus set the number of outputs that you want, using \"Number of
contrasts\". To convert a single EV into a Z statistic image, set it's
contrast value to 1 and all others to 0. Thus the simplest design,
with one EV only, has just one contrast vector, and only one entry in
this contrast vector; 1. To add more contrast vectors, increase the
\"Number of contrasts\". To compare two EVs, for example, to subtract
one stimulus type (EV1) from another type (EV2), set EV1's contrast
value to -1 and EV2's to 1. A Z statistic image will be generated
according to this request."

$w.bhelp bind $fmri(contrastsf).num.ftests -msg "F-tests enable you to investigate several contrasts at the same time,
for example to see whether any of them (or any combination of them) is
significantly non-zero. Also, the F-test allows you to compare the
contribution of each contrast to the model and decide on significant
and non-significant ones.

One example of F-test usage is if a particular stimulation is to be
represented by several EVs, each with the same input function
(e.g. square wave or custom timing) but all with different HRF
convolutions - i.e. several \"basis functions\". Putting all relevant
resulting parameter estimates together into an F-test allows the
complete fit to be tested against zero without having to specify the
relative weights of the basis functions (as one would need to do with
a single contrast). So - if you had three basis functions (EVs 1,2 and
3) the wrong way of combining them is a single (T-test) contrast of
\[1 1 1\]. The right way is to make three contrasts \[1 0 0\] \[0 1 0\] and
\[0 0 1\] and enter all three contrasts into an F-test. As
described above, FEAT will automatically do this for you if you set up
contrasts for \"original EVs\" instead of \"real EVs\".

You can carry out as many F-tests as you like. Each test includes the
particular contrasts that you specify by clicking on the appropriate
buttons."

#}}}
    #{{{ setup EVs

if { $fmri(analysis) != 4 } {

    set fmri(evsf) [ $NB subwidget evs ]

    set ev_description ""
    if { $fmri(level) == 1 } {
	set ev_description " original"
    }
    tixControl $fmri(evsf).evs -label " Number of$ev_description EVs " -variable fmri(evs_orig) -step 1 -min 1 -command "feat5:setup_model_update_evs $w $fmri(evsf) 1"
    $w.bhelp bind $fmri(evsf).evs -msg "The basic number of explanatory variables in the design matrix; this
means the number of different effects that you wish to model - one for
each modelled stimulation type, and one for each modelled confound.

For first-level analyses, it is common for the final design matrix to
have a greater number of \"real EVs\" than this \"original\" number; for
example, when using basis functions, each \"original EV\" gives rise
to several \"real EVs\"."
    pack $fmri(evsf).evs -in $fmri(evsf) -padx 2 -pady 2 -side top -anchor w

    if { $fmri(level) == 1 } {
	tixNoteBook $fmri(evsf).evsnb
	pack $fmri(evsf).evsnb -in $fmri(evsf) -padx 2 -pady 2 -side top -anchor w
    }

    feat5:setup_model_update_evs $w $fmri(evsf) 0 1

} else {
    feat5:setup_model_update_evs $w 0 0 0
}

#}}}
    #{{{ setup buttons

frame $w0.btns
pack $w0.btns -padx 2 -pady 2 -in $w0 -side bottom

frame $w0.btns.b -relief raised -borderwidth 2
pack $w0.btns.b -in $w0.btns -side bottom -fill x -padx 2 -pady 2

button $w0.btns.b.view -command "feat5:setup_model_preview $w" -text "View design"
$w.bhelp bind $w0.btns.b.view -msg $fmri(design_help)

set fmri(cov_help) "This is a graphical representation of the covariance of the design matrix and the efficiency of the design/contrasts. Of most practical importance are the values in the lower part of the window, showing the estimability of the contrasts.


The first matrix shows the absolute value of the normalised correlation of each EV with each EV. If a design is well-conditioned (i.e. not approaching rank deficiency) then the diagonal elements should be white and all others darker.

So - if there are any very bright elements off the diagonal, you can immediately tell which EVs are too similar to each other - for example, if element \[1,3\] (and \[3,1\]) is bright then columns 1 and 3 in the design matrix are possibly too similar.

Note that this includes all real EVs, including any added temporal derivatives, basis functions, etc.

The second matrix shows a similar thing after the design matrix has been run through SVD (singular value decomposition). All non-diagonal elements will be zero and the diagonal elements are given by the eigenvalues of the SVD, so that a poorly-conditioned design is obvious if any of the diagonal elements are black.


In the lower part of the window, for each requested contrast, that contrast's efficiency/estimability is shown. This is formulated as the strength of the signal required in order to detect a statistically significant result for this contrast. For example, in FMRI data and with a single regressor, this shows the BOLD % signal change required. In the case of a differential contrast, it shows the required difference in BOLD signal between conditions.

This \"Effect Required\" depends on the design matrix, the contrast values, the statistical significance level chosen, and the noise level in the data (see the \"Misc\" tab in the main FEAT GUI). The lower the effect required, the more easily estimable is a contrast, i.e. the more efficient is the design.

Note that this does not tell you everything that there is to know about paradigm optimisation. For example, all things being equal, event-related designs tend to give a smaller BOLD effect than block designs - the efficiency estimation made here cannot take that kind of effect into account!"

button $w0.btns.b.acview -command "feat5:setup_model_acpreview $w" -text "Efficiency"
$w.bhelp bind $w0.btns.b.acview -msg $fmri(cov_help)

button $w0.btns.b.cancel -command "feat5:setup_model_destroy $w $w0" -text "Done"
pack $w0.btns.b.view $w0.btns.b.acview -in $w0.btns.b -side left -expand yes -padx 3 -pady 3 -fill y
if { $fmri(infeat) } {
    pack $w0.btns.b.cancel -in $w0.btns.b -side left -expand yes -padx 3 -pady 3 -fill y
}

#}}}
}

#{{{ feat5:setup_model_vars_simple

proc feat5:setup_model_vars_simple { w } {
    global fmri

    set fmri(filmsetup) 1

    if { $fmri(level) == 1 } {

	set fmri(con_mode) orig
	set fmri(con_mode_old) orig

	if { $fmri(wizard_type) == 1 } {
	    #{{{ setup wizard type 1

	set fmri(evs_orig) 1
	set fmri(evs_real) 2

	set fmri(evtitle1) ""
	set fmri(shape1) 0
	set fmri(skip1) 0
	set fmri(off1) $fmri(r_count)
	set fmri(on1) $fmri(a_count)
	set fmri(phase1) 0
	set fmri(stop1) -1

	set fmri(convolve1) $fmri(default_convolve)
	set fmri(convolve_phase1) $fmri(default_convolve_phase)
	set fmri(gammasigma1) $fmri(default_gammasigma)
	set fmri(gammadelay1) $fmri(default_gammadelay)
	set fmri(tempfilt_yn1) 1
	set fmri(deriv_yn1) $fmri(default_deriv_yn)

	set fmri(ortho1.0) 0
	set fmri(ortho1.1) 0

	set fmri(ncon_orig) 1
	set fmri(con_orig1.1) 1
	set fmri(conpic_orig.1) 1
	set fmri(conname_orig.1) ""

	set fmri(nftests_orig) 0

	set fmri(paradigm_hp) [ expr $fmri(r_count) + $fmri(a_count) ]

#}}}
	} elseif { $fmri(wizard_type) == 2 } {
	    #{{{ setup wizard type 2

	set fmri(evs_orig) 2
	set fmri(evs_real) 4

	set fmri(evtitle1) "A"
	set fmri(evtitle2) "B"
	set fmri(shape1) 0
	set fmri(shape2) 0
	set fmri(skip1) 0
	set fmri(skip2) 0
	set fmri(off1) [ expr $fmri(r_count) * 2 + $fmri(b_count) ]
	set fmri(on1) $fmri(a_count)
	set fmri(off2) [ expr $fmri(r_count) * 2 + $fmri(a_count) ]
	set fmri(on2) $fmri(b_count)
	set fmri(phase1) [ expr $fmri(r_count) + $fmri(b_count) ]
	set fmri(phase2) 0
	set fmri(stop1) -1
	set fmri(stop2) -1
	set fmri(convolve1) $fmri(default_convolve)
	set fmri(convolve2) $fmri(default_convolve)
	set fmri(convolve_phase1) $fmri(default_convolve_phase)
	set fmri(convolve_phase2) $fmri(default_convolve_phase)
	set fmri(gammasigma1) $fmri(default_gammasigma)
	set fmri(gammasigma2) $fmri(default_gammasigma)
	set fmri(gammadelay1) $fmri(default_gammadelay)
	set fmri(gammadelay2) $fmri(default_gammadelay)
	set fmri(tempfilt_yn1) 1
	set fmri(tempfilt_yn2) 1
	set fmri(deriv_yn1) $fmri(default_deriv_yn)
	set fmri(deriv_yn2) $fmri(default_deriv_yn)

	set fmri(ortho1.0) 0
	set fmri(ortho1.1) 0
	set fmri(ortho1.2) 0
	set fmri(ortho2.0) 0
	set fmri(ortho2.1) 0
	set fmri(ortho2.2) 0
	
	set fmri(ncon_orig) 4
	set fmri(con_orig1.1) 1
	set fmri(con_orig1.2) 0
	set fmri(con_orig2.1) 0
	set fmri(con_orig2.2) 1
	set fmri(con_orig3.1) 1
	set fmri(con_orig3.2) -1
	set fmri(con_orig4.1) -1
	set fmri(con_orig4.2) 1

	set fmri(conpic_orig.1) 1
	set fmri(conpic_orig.2) 1
	set fmri(conpic_orig.3) 1
	set fmri(conpic_orig.4) 1

	set fmri(conname_orig.1) "A"
	set fmri(conname_orig.2) "B"
	set fmri(conname_orig.3) "A>B"
	set fmri(conname_orig.4) "B>A"

	set fmri(nftests_orig) 1
	set fmri(ftest_orig1.1) 1
	set fmri(ftest_orig1.2) 1
	set fmri(ftest_orig1.3) 0
	set fmri(ftest_orig1.4) 0

	set fmri(paradigm_hp) [ expr  ( 2 * $fmri(r_count) ) + $fmri(a_count) + $fmri(b_count) ]

#}}}
	} else {
	    #{{{ setup wizard type 3

	set fmri(evs_orig) 3
	set fmri(evs_real) 3

	set fmri(evtitle1) "c-t"
	set fmri(evtitle2) "BOLD"
	set fmri(evtitle3) "c-t act"
	set fmri(shape1) 0
	set fmri(shape2) 0
	set fmri(shape3) 4
	set fmri(skip1) 0
	set fmri(skip2) 0
	set fmri(off1) $fmri(tr)
	set fmri(on1) $fmri(tr)
	set fmri(off2) $fmri(r_count)
	set fmri(on2) $fmri(a_count)
	set fmri(phase1) 0
	set fmri(phase2) 0
	set fmri(stop1) -1
	set fmri(stop2) -1
	set fmri(convolve1) 0
	set fmri(convolve2) $fmri(default_convolve)
	set fmri(convolve3) 0
	set fmri(convolve_phase2) $fmri(default_convolve_phase)
	set fmri(convolve_phase3) 0
	set fmri(gammasigma2) $fmri(default_gammasigma)
	set fmri(gammadelay2) $fmri(default_gammadelay)

        set fmri(interactions3.1) 1
        set fmri(interactionsd3.1) 1
        set fmri(interactions3.2) 1
        set fmri(interactionsd3.2) 0

	set fmri(tempfilt_yn1) 1
	set fmri(tempfilt_yn2) 1
	set fmri(tempfilt_yn3) 1
	set fmri(deriv_yn1) 0
	set fmri(deriv_yn2) 0
	set fmri(deriv_yn3) 0

	set fmri(ortho1.0) 0
	set fmri(ortho1.1) 0
	set fmri(ortho1.2) 0
	set fmri(ortho1.3) 0
	set fmri(ortho2.0) 0
	set fmri(ortho2.1) 0
	set fmri(ortho2.2) 0
	set fmri(ortho2.3) 0
	set fmri(ortho3.0) 0
	set fmri(ortho3.1) 0
	set fmri(ortho3.2) 0
	set fmri(ortho3.3) 0
	
	set fmri(ncon_orig) 6
	set fmri(con_orig1.1) 0
	set fmri(con_orig1.2) 0
	set fmri(con_orig1.3) 1
	set fmri(con_orig2.1) 0
	set fmri(con_orig2.2) 0
	set fmri(con_orig2.3) -1
	set fmri(con_orig3.1) 0
	set fmri(con_orig3.2) 1
	set fmri(con_orig3.3) 0
	set fmri(con_orig4.1) 0
	set fmri(con_orig4.2) -1
	set fmri(con_orig4.3) 0
	set fmri(con_orig5.1) 1
	set fmri(con_orig5.2) 0
	set fmri(con_orig5.3) 0
	set fmri(con_orig6.1) -1
	set fmri(con_orig6.2) 0
	set fmri(con_orig6.3) 0

	set fmri(conpic_orig.1) 1
	set fmri(conpic_orig.2) 1
	set fmri(conpic_orig.3) 1
	set fmri(conpic_orig.4) 1
	set fmri(conpic_orig.5) 1
	set fmri(conpic_orig.6) 1

	set fmri(conname_orig.1) "perfusion activation"
	set fmri(conname_orig.2) "-perfusion activation"
	set fmri(conname_orig.3) "BOLD"
	set fmri(conname_orig.4) "-BOLD"
	set fmri(conname_orig.5) "control-tag baseline"
	set fmri(conname_orig.6) "control-tag -baseline"

	set fmri(nftests_orig) 0

	set fmri(paradigm_hp) [ expr $fmri(r_count) + $fmri(a_count) ]

#}}}
	}

    } else {

	set fmri(tr) 3

	if { $fmri(wizard_type) == 1 } {
	    #{{{ setup higher-level starting values

set fmri(evs_orig) 1
set fmri(evs_real) 1

set fmri(custom1) dummy
set fmri(shape1) 2
set fmri(convolve1) 0    
set fmri(tempfilt_yn1) 0
set fmri(deriv_yn1) 0
set fmri(ortho1.0) 0
for { set i 1 } { $i <= $fmri(multiple) } { incr i 1 } {
    set fmri(evg${i}.1) 1
    set fmri(groupmem.${i}) 1
}

set fmri(ncon_real) 1
set fmri(con_mode) real
set fmri(con_mode_old) real
set fmri(con_real1.1) 1
set fmri(conpic_real.1) 1
set fmri(nftests_real) 0
set fmri(conname_real.1) "group mean"

#}}}
	} elseif { $fmri(wizard_type) == 2 } {
	    #{{{ setup higher-level starting values

set fmri(evs_orig) 2
set fmri(evs_real) 2

set fmri(evtitle1) "group A"
set fmri(evtitle2) "group B"
set fmri(custom1) dummy
set fmri(custom2) dummy
set fmri(shape1) 2
set fmri(shape2) 2
set fmri(convolve1) 0    
set fmri(convolve2) 0    
set fmri(convolve_phase1) 0    
set fmri(convolve_phase2) 0    
set fmri(tempfilt_yn1) 0
set fmri(tempfilt_yn2) 0
set fmri(deriv_yn1) 0
set fmri(deriv_yn2) 0
set fmri(ortho1.0) 0
set fmri(ortho1.1) 0
set fmri(ortho1.2) 0
set fmri(ortho2.0) 0
set fmri(ortho2.1) 0
set fmri(ortho2.2) 0
for { set i 1 } { $i <=  $fmri(multiple) } { incr i 1 } {
    if { $i <= $fmri(a_count) } {
	set fmri(evg${i}.1) 1
	set fmri(evg${i}.2) 0
	set fmri(groupmem.${i}) 1
    } else {
	set fmri(evg${i}.1) 0
	set fmri(evg${i}.2) 1
	set fmri(groupmem.${i}) 2
    }
}

set fmri(ncon_real) 4
set fmri(con_mode) real
set fmri(con_mode_old) real
set fmri(conpic_real.1) 1
set fmri(conpic_real.2) 1
set fmri(conpic_real.3) 1
set fmri(conpic_real.4) 1
set fmri(con_real1.1) 1
set fmri(con_real1.2) -1
set fmri(con_real2.1) -1
set fmri(con_real2.2) 1
set fmri(con_real3.1) 1
set fmri(con_real3.2) 0
set fmri(con_real4.1) 0
set fmri(con_real4.2) 1
set fmri(conname_real.1) "group A > group B"
set fmri(conname_real.2) "group B > group A"
set fmri(conname_real.3) "group A mean"
set fmri(conname_real.4) "group B mean"
set fmri(nftests_real) 0

#}}}
	} else {
	    #{{{ setup higher-level starting values

set nsubjects [ expr $fmri(multiple) / 2 ]

set fmri(evs_orig) [ expr $nsubjects + 1 ]
set fmri(evs_real) $fmri(evs_orig)

for { set i 1 } { $i <=  $fmri(evs_orig) } { incr i 1 } {
    set fmri(custom${i}) dummy
    set fmri(shape${i}) 2
    set fmri(convolve${i}) 0    
    set fmri(convolve_phase${i}) 0    
    set fmri(tempfilt_yn${i}) 0
    set fmri(deriv_yn${i}) 0
    set fmri(evtitle${i}) "s[ expr $i - 1 ]"
    for { set j 0 } { $j <=  $fmri(evs_orig) } { incr j 1 } {
	set fmri(ortho${i}.${j}) 0
    }
}

set fmri(evtitle1) "A > B"

for { set i 1 } { $i <=  $fmri(multiple) } { incr i 1 } {
    set fmri(groupmem.${i}) 1

    if { $i <= $nsubjects } {
	set fmri(evg${i}.1) 1
    } else {
	set fmri(evg${i}.1) -1
    }

    for { set j 2 } { $j <= $fmri(evs_orig) } { incr j 1 } {
	set fmri(evg${i}.${j}) 0	
    }
}

for { set i 1 } { $i <= $nsubjects } { incr i 1 } {
    set fmri(evg${i}.[ expr 1 + $i ]) 1
    set fmri(evg[ expr $i + $nsubjects ].[ expr 1 + $i ]) 1
}

set fmri(ncon_real) 2
set fmri(con_mode) real
set fmri(con_mode_old) real
set fmri(conpic_real.1) 1
set fmri(conpic_real.2) 1

set fmri(con_real1.1) 1
set fmri(con_real2.1) -1
for { set i 2 } { $i <= $fmri(evs_orig) } { incr i 1 } {
    set fmri(con_real1.${i}) 0
    set fmri(con_real2.${i}) 0
}

set fmri(conname_real.1) "condition A > B"
set fmri(conname_real.2) "condition B > A"
set fmri(nftests_real) 0

#}}}
	}
    }
}

#}}}
#{{{ feat5:setup_model_update_ev_i

proc feat5:setup_model_update_ev_i { w w0 i force dummy } {
    global fmri

    if { $fmri(level)==1 } {

	#{{{ basic shape/timings stuff

pack forget $w0.evsnb.skip$i $w0.evsnb.off$i $w0.evsnb.on$i $w0.evsnb.phase$i $w0.evsnb.stop$i $w0.evsnb.period$i $w0.evsnb.nharmonics$i $w0.evsnb.custom$i $w0.evsnb.interaction$i

if { $fmri(shape$i) == 0 } {

    pack $w0.evsnb.skip$i $w0.evsnb.off$i $w0.evsnb.on$i $w0.evsnb.phase$i $w0.evsnb.stop$i -in $w0.evsnb.timings$i -padx 5 -pady 2 -side top -anchor w

} elseif { $fmri(shape$i) == 1 } {

    pack $w0.evsnb.skip$i $w0.evsnb.period$i $w0.evsnb.phase$i $w0.evsnb.nharmonics$i $w0.evsnb.stop$i -in $w0.evsnb.timings$i -padx 5 -pady 2 -side top -anchor w

} elseif { $fmri(shape$i) == 2 || $fmri(shape$i) == 3 } {

    pack $w0.evsnb.custom$i -in $w0.evsnb.timings$i -padx 5 -pady 2 -side top -anchor w

} elseif { $fmri(shape$i) == 4 } {

    pack $w0.evsnb.interaction$i -in $w0.evsnb.timings$i -padx 5 -pady 2 -side top -anchor w
    set selected 0
    for { set j 1 } { $j < $i } { incr j 1 } {
	set selected [ expr $selected + $fmri(interactions${i}.$j) ]
    }
    if { $selected < 2 } {
	if { !$fmri(interactions${i}.1) } {
	    set fmri(interactions${i}.1) 1
	} else {
	    set fmri(interactions${i}.2) 1
	}
    }

}

if { $force } {

    if { $fmri(shape$i) == 1 } {
	set fmri(phase$i) -6
    } else {
	set fmri(phase$i) 0
    }

    if { $fmri(shape$i) != 4 } {
	set fmri(convolve$i) $fmri(default_convolve)
	set fmri(tempfilt_yn$i) 1
    }
}

#}}}

	if { $fmri(shape$i) != 1 && $fmri(shape$i) != 4 && $fmri(shape$i) != 10 } {
	    pack $w0.evsnb.conv$i -in $fmri(modelf$i) -after $w0.evsnb.timings$i -padx 5 -pady 2 -side top -anchor w
	    pack forget $w0.evsnb.convolve_phase$i $w0.evsnb.gausssigma$i $w0.evsnb.gaussdelay$i \
		    $w0.evsnb.gammasigma$i $w0.evsnb.gammadelay$i $w0.evsnb.bfcustom$i $w0.evsnb.bfcustomlabel$i $w0.evsnb.basisfnum$i $w0.evsnb.basisfwidth$i 
	    if { $fmri(convolve$i) > 0 } {
		pack $w0.evsnb.convolve_phase$i -in $w0.evsnb.conv$i -padx 5 -pady 2 -side top -anchor w
		if { $fmri(convolve$i) == 1 } {
		    pack $w0.evsnb.gausssigma$i $w0.evsnb.gaussdelay$i -in $w0.evsnb.conv$i -padx 5 -pady 2 -side top -anchor w
		} elseif { $fmri(convolve$i) == 2 } {
		    pack $w0.evsnb.gammasigma$i $w0.evsnb.gammadelay$i -in $w0.evsnb.conv$i -padx 5 -pady 2 -side top -anchor w
		} elseif { $fmri(convolve$i) > 3 && $fmri(convolve$i) < 7 } {
		    pack $w0.evsnb.basisfnum$i $w0.evsnb.basisfwidth$i -in $w0.evsnb.conv$i -padx 5 -pady 2 -side top -anchor w
		} elseif { $fmri(convolve$i) == 7 } {
		    pack $w0.evsnb.bfcustom$i $w0.evsnb.bfcustomlabel$i -in $w0.evsnb.conv$i -padx 5 -pady 2 -side top -anchor w
		}
	    }
	} else {
	    pack forget $w0.evsnb.conv$i
	    set fmri(convolve$i) 0
	}

	#{{{ tempfilt

if { $fmri(shape$i) != 10 } {
    pack $w0.evsnb.tempfilt$i -in $fmri(modelf$i) -padx 5 -pady 2 -side bottom -anchor w
} else { 
    pack forget $w0.evsnb.tempfilt$i
}

#}}}

	#{{{ orthogonalise

if { $fmri(evs_orig) > 1 && $fmri(shape$i) != 10 } {

    if { ! [ winfo exists $w0.evsnb.ortho$i ] } {

	frame $w0.evsnb.ortho$i

	checkbutton $w0.evsnb.ortho$i.0 -variable fmri(ortho${i}.0) \
		-command "feat5:setup_model_update_ev_i $w $w0 $i 0 0"

	pack $w0.evsnb.ortho$i.0 -in $w0.evsnb.ortho$i -side left

	$w.bhelp bind $w0.evsnb.ortho$i -msg "Orthogonalising an EV with respect to other EVs means that it is
completely independent of the other EVs, i.e. contains no component
related to them. Most sensible designs are already in this form - all
EVs are at least close to being orthogonal to all others. However,
this may not be the case; you can use this facility to force an EV to
be orthogonal to some or all other EVs. This is achieved by
subtracting from the EV that part which is related to the other EVs
selected here.

An example use would be if you had another EV which was a
constant height spike train, and the current EV is derived from this
other one, but with a linear increase in spike height imposed, to
model an increase in response during the experiment for any
reason. You would not want the current EV to contain any component of
the constant height EV, so you would orthogonalise the current EV wrt
the other."

        pack $w0.evsnb.ortho$i -in $fmri(modelf$i) -after $w0.evsnb.tempfilt$i -padx 5 -pady 2 -side top -anchor w
    }

    for { set j 1 } { $j > 0 } { incr j 1 } {
	if { [ winfo exists $w0.evsnb.ortho$i.$j ] } {
	    destroy $w0.evsnb.ortho$i.$j
	} elseif { $j != $i } {
	    set j -10
	}
    }

    if { $fmri(ortho${i}.0) == 1 } {
	$w0.evsnb.ortho$i.0 configure -text "Orthogonalise    wrt EVs "
	for { set j 1 } { $j <= $fmri(evs_orig) } { incr j 1 } {
	    if { $j != $i } {
		checkbutton $w0.evsnb.ortho$i.$j -text "$j " -variable fmri(ortho${i}.$j)
		pack $w0.evsnb.ortho$i.$j -in $w0.evsnb.ortho$i -side left -padx 0
	    }
	}
    } else {
	$w0.evsnb.ortho$i.0 configure -text "Orthogonalise"
    }

} else {

    if { [ winfo exists $w0.evsnb.ortho$i ] } {
	destroy $w0.evsnb.ortho$i $w0.evsnb.ortho$i.0
	set fmri(ortho${i}.0) 0
    }
}

#}}}

	#{{{ deriv

if { ( $fmri(shape$i) != 1 && $fmri(convolve$i) < 4 ) || $fmri(shape$i) == 10 } {
    pack $w0.evsnb.deriv$i -in $fmri(modelf$i) -padx 5 -pady 2 -side bottom -anchor w
} else { 
    set fmri(deriv_yn$i) 0
    pack forget $w0.evsnb.deriv$i
}

#}}}
    }

    feat5:setup_model_update_contrasts $w 0
}

#}}}
#{{{ feat5:setup_model_update_evs

proc feat5:setup_model_update_evs { w w0 dummy update_gui } {
    global fmri PWD

    #{{{ initialise variables

    for { set i 1 } { $i <= $fmri(evs_orig) } { incr i 1 } {

	for { set j 0 } { $j <= $fmri(evs_orig) } { incr j 1 } {
	    if { ! [ info exists fmri(ortho${i}.$j) ] } {
		set fmri(ortho${i}.$j) 0
	    }
	}

	if { $fmri(level) == 1 } {
	    #{{{ initialise variables

if { ! [ info exists fmri(evtitle$i) ] } {
    set fmri(evtitle$i) ""
}

if { ! [ info exists fmri(shape$i) ] } {
    set fmri(shape$i) 0
}

if { ! [ info exists fmri(skip$i) ] } {
    set fmri(skip$i) 0
}

if { ! [ info exists fmri(off$i) ] } {
    set fmri(off$i) 30
}

if { ! [ info exists fmri(on$i) ] } {
    set fmri(on$i) 30
}

if { ! [ info exists fmri(phase$i) ] } {
    set fmri(phase$i) 0
}

if { ! [ info exists fmri(stop$i) ] } {
    set fmri(stop$i) -1
}

if { ! [ info exists fmri(period$i) ] } {
    set fmri(period$i) 60
}

if { ! [ info exists fmri(nharmonics$i) ] } {
    set fmri(nharmonics$i) 0
}

if { ! [ info exists fmri(convolve$i) ] } {
    set fmri(convolve$i) $fmri(default_convolve)
}

if { ! [ info exists fmri(convolve_phase$i) ] } {
    set fmri(convolve_phase$i) $fmri(default_convolve_phase)
}

if { ! [ info exists fmri(gausssigma$i) ] } {
    set fmri(gausssigma$i) $fmri(default_gausssigma)
}

if { ! [ info exists fmri(gaussdelay$i) ] } {
    set fmri(gaussdelay$i) $fmri(default_gaussdelay)
}

if { ! [ info exists fmri(gammasigma$i) ] } {
    set fmri(gammasigma$i) $fmri(default_gammasigma)
}

if { ! [ info exists fmri(gammadelay$i) ] } {
    set fmri(gammadelay$i) $fmri(default_gammadelay)
}

if { ! [ info exists fmri(bfcustom$i) ] } {
    set fmri(bfcustom$i) $fmri(default_bfcustom) 
}

if { ! [ info exists fmri(basisfnum$i) ] } {
    set fmri(basisfnum$i) 3
}

if { ! [ info exists fmri(basisfwidth$i) ] } {
    set fmri(basisfwidth$i) 15
}

if { ! [ info exists fmri(tempfilt_yn$i) ] } {
    set fmri(tempfilt_yn$i) 1
}

if { ! [ info exists fmri(deriv_yn$i) ] } {
    set fmri(deriv_yn$i) $fmri(default_deriv_yn)
}

if { ! [ info exists fmri(interactions${i}.1) ] && $i > 2 } {
    for { set j 1 } { $j < $i } { incr j 1 } {
	if { $j < 3 } {
	    set fmri(interactions${i}.$j) 1
	} else {
	    set fmri(interactions${i}.$j) 0
	}
	set fmri(interactionsd${i}.$j) 0
    }
}

#}}}
	} else {
	    set fmri(custom$i) dummy
	    set fmri(shape$i) 2
	    set fmri(convolve$i) 0    
	    set fmri(convolve_phase$i) 0
	    set fmri(tempfilt_yn$i) 0
	    set fmri(deriv_yn$i) 0
	    for { set j 1 } { $j <= $fmri(evs_orig) } { incr j 1 } {
		if { ! [ info exists fmri(evg$i.$j) ] } {
		    set fmri(evg$i.$j) 0
		}
	    }
	}
    }

#}}}

    if { $update_gui } {

	if { $fmri(level) == 1 } {
	    #{{{ setup EV notebook for level=1

	set LWIDTH 15

	for { set i 1 } { $i > 0 } { incr i 1 } {
	    if { $i <= $fmri(evs_orig) } {

		if { ! [ info exists fmri(modelf$i) ] || ! [ winfo exists $fmri(modelf$i) ] } {
		    #{{{ setup EV tab i

$w0.evsnb add ev$i -label "$i"
set fmri(modelf$i) [ $w0.evsnb subwidget ev$i ]

if { $fmri(level) == 1 } {

    #{{{ EV name

frame $w0.evsnb.evtitle$i

label $w0.evsnb.evtitle${i}.label -text "EV name "

entry $w0.evsnb.evtitle${i}.entry -textvariable fmri(evtitle$i) -width 7

$w.bhelp bind $w0.evsnb.evtitle$i -msg "If you wish, enter a title for EV $i here."

pack $w0.evsnb.evtitle${i}.label $w0.evsnb.evtitle${i}.entry -in $w0.evsnb.evtitle$i -side left

#}}}

    #{{{ basic shape

    set grot $fmri(shape$i)

    tixOptionMenu $w0.evsnb.shape$i -label "Basic shape: " -variable fmri(shape$i)
    $w0.evsnb.shape$i add command 10 -label "Empty (all zeros)"
    $w0.evsnb.shape$i add command 0 -label "Square"
    $w0.evsnb.shape$i add command 1 -label "Sinusoid"
    $w0.evsnb.shape$i add command 2 -label "Custom (1 entry per volume)"
    $w0.evsnb.shape$i add command 3 -label "Custom (3 column format)"
    if { $i > 2 } {
	$w0.evsnb.shape$i add command 4 -label "Interaction"
    }

    set fmri(shape$i) $grot

$w.bhelp bind $w0.evsnb.shape$i -msg "Choose the basic shape of the waveform that describes the stimulus or confound that you wish to model. The basic waveform should be exactly in time with the applied stimulation, i.e., not lagged at all. This is because the measured (time-series) response will be delayed with respect to the stimulation, and this delay is modelled in the design matrix by convolution of the basic waveform with a suitable haemodynamic response function (see appropriate bubble-help).

If you need an EV to be ignored, choose \"Empty (all zeros)\". You are most likely to want to do this if you want the EVs to all have the same meaning for multiple runs, but in some runs one or more EVs contain no events of the relevant type. Note that in this case you will get a warning about the matrix being rank deficient.

For an on/off (or a regularly-spaced single-event) experiment choose a \"square\" wave. To model single-event experiments with this method, the \"On\" periods will probably be small - e.g., 1s or even less.

For sinusoidal modelling choose the \"Sinusoid\" option and select the number of \"Harmonics\" (or overtones) that you want to add to the fundamental frequency.

For a single-event experiment with irregular timing for the stimulations, a custom file can be used. With \"Custom (1 entry per volume)\", you specify a single value for each timepoint. The custom file should be a raw text file, and should be a list of numbers, separated by spaces or newlines, with one number for each volume (after subtracting the number of deleted images). These numbers can either all be 0s and 1s, or can take a range of values. The former case would be appropriate if the same stimulus was applied at varying time points; the latter would be appropriate, for example, if recorded subject responses are to be inserted as an effect to be modelled. Note that it may or may not be appropriate to convolve this particular waveform with an HRF - in the case of single-event, it is.

For even finer control over the input waveform, choose \"Custom (3 column format)\". In this case the custom file consists of triplets of numbers; you can have any number of triplets. Each triplet describes a short period of time and the value of the model during that time. The first number in each triplet is the onset (in seconds) of the period, the second number is the duration (in seconds) of the period, and the third number is the value of the input during that period. The same comments as above apply, about whether these numbers are 0s and 1s, or vary continuously. The start of the first non-deleted volume correpsonds to t=0.

Note that whilst ALL columns are demeaned before model fitting, neither custom format will get rescaled - it is up to you to make sure that relative scaling between different EVs is sensible. If you double the scaling of values in an EV you will halve the resulting parameter estimate, which will change contrasts of this EV against others.

If you select \"Interaction\" then the current EV is modelled as an interaction between other EVs, and is normally used to create a third EV from two existing EVs, to model the nonlinear interaction between two different conditions. On the line of buttons marked \"Between EVs\" you select which other EVs to interact to form the current one. The selected EVs then get multiplied together to form the current EV. Normally they are multiplied after (temporarily) shifting their values so that the minimum of each EV is zero; however, if you click on the relevant \"Demean EV\" button, they are zero-centred instead (this demeaning is normally only used when forming perfusion FMRI models)."

#}}}
    #{{{ timings / custom name

frame $w0.evsnb.timings$i

tixControl $w0.evsnb.skip$i -variable fmri(skip$i) -step 1 -min 0 -selectmode immediate -label "    Skip (s)"
$w0.evsnb.skip$i subwidget label config -width $LWIDTH
$w.bhelp bind $w0.evsnb.skip$i -msg "The initial period (seconds) before the waveform commences."

tixControl $w0.evsnb.off$i -variable fmri(off$i) -step 1 -min 0 -selectmode immediate -label "    Off (s)"
$w0.evsnb.off$i subwidget label config -width $LWIDTH
$w.bhelp bind $w0.evsnb.off$i -msg "The duration (seconds) of the \"Off\" periods in the square wave."

tixControl $w0.evsnb.on$i -variable fmri(on$i) -step 1 -min 0 -selectmode immediate -label "    On (s)"
$w0.evsnb.on$i subwidget label config -width $LWIDTH
$w.bhelp bind $w0.evsnb.on$i -msg "The duration (seconds) of the \"On\" periods in the square wave."

tixControl $w0.evsnb.phase$i -variable fmri(phase$i) -step 1 -selectmode immediate -label "    Phase (s)"
$w0.evsnb.phase$i subwidget label config -width $LWIDTH
$w.bhelp bind $w0.evsnb.phase$i -msg "The phase shift (seconds) of the waveform. By default, after the \"Skip\" period, the square wave\nstarts with a full \"Off\" period and the \"Sinusoid\" starts by falling from zero. However, the\nwave can be brought forward in time according to the phase shift."

tixControl $w0.evsnb.stop$i -variable fmri(stop$i) -step 1 -min -1 -selectmode immediate -label "    Stop after (s)"
$w0.evsnb.stop$i subwidget label config -width $LWIDTH
$w.bhelp bind $w0.evsnb.stop$i -msg "The active duration (seconds) of the waveform, starting\nafter the \"Skip\" period. \"-1\" means do not stop."

tixControl $w0.evsnb.period$i -variable fmri(period$i) -step 1 -min 0 -selectmode immediate -label "    Period (s)"
$w0.evsnb.period$i subwidget label config -width $LWIDTH
$w.bhelp bind $w0.evsnb.period$i -msg "The period (seconds) of the \"Sinusoid\" waveform."

tixControl $w0.evsnb.nharmonics$i -variable fmri(nharmonics$i) -step 1 -min 0 -selectmode immediate -label "    Harmonics" \
	-command "feat5:setup_model_update_contrasts $w"
$w0.evsnb.nharmonics$i subwidget label config -width $LWIDTH
$w.bhelp bind $w0.evsnb.nharmonics$i -msg "How many harmonics (sine waves with periods of half the primary sine
wave, then quarter, etc) would you like?"

FSLFileEntry $w0.evsnb.custom$i \
	-variable fmri(custom$i) \
	-pattern "*" \
	-directory $PWD \
	-label "    Filename" \
	-title "Select an event file" \
	-width 30 \
	-filterhist VARS(history)

if { $i > 2 } {
    frame $w0.evsnb.interaction$i
    label $w0.evsnb.interaction$i.label -text "Between EVs "
    label $w0.evsnb.interaction$i.labeld -text "Demean EV "
    grid $w0.evsnb.interaction$i.label -in $w0.evsnb.interaction$i -column 0 -row 0
    grid $w0.evsnb.interaction$i.labeld -in $w0.evsnb.interaction$i -column 0 -row 1
    for { set j 1 } { $j < $i } { incr j 1 } {
	checkbutton $w0.evsnb.interaction$i.$j -variable fmri(interactions${i}.$j) -text "$j " -command "feat5:setup_model_update_ev_i $w $w0 $i 0 0"
	grid $w0.evsnb.interaction$i.$j -in $w0.evsnb.interaction$i -column $j -row 0
	checkbutton $w0.evsnb.interaction$i.d$j -variable fmri(interactionsd${i}.$j) -text "$j "
	grid $w0.evsnb.interaction$i.d$j -in $w0.evsnb.interaction$i -column $j -row 1
    }
}

#}}}

    #{{{ convolution

frame $w0.evsnb.conv$i

set tmpval $fmri(convolve$i)
tixOptionMenu $w0.evsnb.convolve$i -label "Convolution: " -variable fmri(convolve$i)
$w0.evsnb.convolve$i add command 0 -label "None"
$w0.evsnb.convolve$i add command 1 -label "Gaussian"
$w0.evsnb.convolve$i add command 2 -label "Gamma"
$w0.evsnb.convolve$i add command 3 -label "Double-Gamma HRF"
$w0.evsnb.convolve$i add command 7 -label "Optimal/custom basis functions"
$w0.evsnb.convolve$i add command 4 -label "Gamma basis functions"
$w0.evsnb.convolve$i add command 5 -label "Sine basis functions"
$w0.evsnb.convolve$i add command 6 -label "FIR basis functions"
set fmri(convolve$i) $tmpval

$w.bhelp bind $w0.evsnb.convolve$i -msg "The form of the HRF (haemodynamic response function) convolution that
will be applied to the basic waveform. This blurs and delays the
original waveform, in an attempt to match the difference between the
input function (original waveform, i.e., stimulus waveform) and the
output function (measured FMRI haemodynamic response). 

If the original waveform is already in an appropriate form, e.g., was
sampled from the data itself, \"None\" should be selected.

The next three options are all somewhat similar blurring and delaying
functions. \"Gaussian\" is simply a Gaussian kernel, whose width and
lag can be altered. \"Gamma\" is a Gamma variate (in fact a
normalisation of the probability density function of the Gamma
function); again, width and lag can be altered. \"Double-Gamma HRF\" is a
preset function which is a mixture of two Gamma functions - a standard
positive function at normal lag, and a small, delayed, inverted Gamma,
which attempts to model the late undershoot. 

The remaining convolution options setup different \"basis
functions\". This means that the original EV waveform will get
convolved by a \"basis set\" of related but different convolution
kernels. By default, an \"original EV\" will generate a set of \"real
EVs\", one for each basis function.

The \"Optimal/custom\" option allows you to use a customised set of
basis functions, setup in a plain text file with one column for each
basis function, sampled at the temporal resolution of 0.05s. The main
point of this option is to allow the use of \"FLOBS\" (FMRIB's Linear
Optimal Basis Set), which is a method for generating a set of basis
functions that has optimal efficiency in covering the range of likely
HRF shapes actually found in your data. You can either use the default
FLOBS set, or use the \"Make_flobs\" GUI on the FEAT \"Utils\" menu to
create your own customised set of FLOBS.

The other basis function options, which will not in general be as good
at fitting the data as FLOBS, are a set of \"Gamma\" variates of
different widths and lags, a set of \"Sine\" waves of differing
frequencies or a set of \"FIR\" (finite-impulse-response) filters
(with FIR the convolution kernel is represented as a set of discrete
fixed-width \"impulses\")."

pack $w0.evsnb.convolve$i -in $w0.evsnb.conv$i -padx 0 -pady 2 -side top -anchor w

#}}}
    #{{{ convolution parameters

tixControl $w0.evsnb.convolve_phase$i -variable fmri(convolve_phase$i) -step 1 -selectmode immediate -label "    Phase (s)"
$w0.evsnb.convolve_phase$i subwidget label config -width $LWIDTH
$w.bhelp bind $w0.evsnb.convolve_phase$i -msg "This sets the \"Phase\" of the convolution - i.e. phase shifts the
convolved time series. Positive values shift the final time series
earlier in time."

tixControl $w0.evsnb.gausssigma$i -variable fmri(gausssigma$i) -step 0.1 -min 0 -selectmode immediate -label "    Sigma (s)"
$w0.evsnb.gausssigma$i subwidget label config -width $LWIDTH
$w.bhelp bind $w0.evsnb.gausssigma$i -msg "This sets the half-width of the Gaussian smoothing of the
input waveform."

tixControl $w0.evsnb.gaussdelay$i -variable fmri(gaussdelay$i) -step 0.1 -min 0 -selectmode immediate -label "    Peak lag (s)"
$w0.evsnb.gaussdelay$i subwidget label config -width $LWIDTH
$w.bhelp bind $w0.evsnb.gaussdelay$i -msg "This sets the peak lag of the Gaussian smoothing of the
input waveform."

tixControl $w0.evsnb.gammasigma$i -variable fmri(gammasigma$i) -step 0.1 -min 0.01 -selectmode immediate -label "    Stddev (s)"
$w0.evsnb.gammasigma$i subwidget label config -width $LWIDTH
$w.bhelp bind $w0.evsnb.gammasigma$i -msg "This sets the half-width of the Gamma smoothing of the
input waveform."

tixControl $w0.evsnb.gammadelay$i -variable fmri(gammadelay$i) -step 0.1 -min 0.01 -selectmode immediate -label "    Mean lag (s)"
$w0.evsnb.gammadelay$i subwidget label config -width $LWIDTH
$w.bhelp bind $w0.evsnb.gammadelay$i -msg "This sets the mean lag of the Gamma smoothing of the
input waveform."

FSLFileEntry $w0.evsnb.bfcustom$i \
	-variable fmri(bfcustom$i) \
	-pattern "*" \
	-directory $PWD \
	-label "    Filename" \
	-title "Select a custom HRF convolution file" \
	-width 30 \
	-filterhist VARS(history) \
        -command "feat5:checkbfcustom $w $i"

label $w0.evsnb.bfcustomlabel$i -text "      (create a custom optimal basis set with Utils->Make_flobs)"

tixControl $w0.evsnb.basisfnum$i -variable fmri(basisfnum$i) -step 1 -min 1 -selectmode immediate \
	-command "feat5:setup_model_update_contrasts $w" -label "    Number"
$w0.evsnb.basisfnum$i subwidget label config -width $LWIDTH
$w.bhelp bind $w0.evsnb.basisfnum$i -msg "This sets the number of basis functions."

tixControl $w0.evsnb.basisfwidth$i -variable fmri(basisfwidth$i) -step 1 -min 1 -selectmode immediate -label "    Window (s)"
$w0.evsnb.basisfwidth$i subwidget label config -width $LWIDTH
$w.bhelp bind $w0.evsnb.basisfwidth$i -msg "This sets the total period over which the basis functions are spread."

#}}}

    #{{{ temporal filtering

checkbutton $w0.evsnb.tempfilt$i -variable fmri(tempfilt_yn$i) -text "Apply temporal filtering" 

$w.bhelp bind $w0.evsnb.tempfilt$i -msg "You should normally apply the same temporal filtering to the model as
you have applied to the data, as the model is designed to look like
the data before temporal filtering was applied. Thus long-time-scale
components in the model will be dealt with correctly."

#}}}
    #{{{ temporal derivative

checkbutton $w0.evsnb.deriv$i -variable fmri(deriv_yn$i) -text "Add temporal derivative" \
	-command "feat5:setup_model_update_contrasts $w 0"

$w.bhelp bind $w0.evsnb.deriv$i -msg "Adding a fraction of the temporal derivative of the blurred original
waveform is equivalent to shifting the waveform slightly in time, in
order to achieve a slightly better fit to the data. Thus adding in the
temporal derivative of a waveform into the design matrix allows a
better fit for the whole model, reducing unexplained noise, and
increasing resulting statistical significances. This option is not
available if you are using basis functions."

#}}}

    pack $w0.evsnb.evtitle$i $w0.evsnb.shape$i $w0.evsnb.timings$i $w0.evsnb.tempfilt$i -in $fmri(modelf$i) -padx 5 -pady 2 -side top -anchor w

    $w0.evsnb.shape$i configure  -command " feat5:setup_model_update_ev_i $w $w0 $i 1"
    $w0.evsnb.convolve$i configure -command " feat5:setup_model_update_ev_i $w $w0 $i 0"

} else {

    for { set j 1 } { $j <= $fmri(multiple) } { incr j 1 } {
	tixControl $w0.evg${i}_${j} -variable fmri(evg${i}.${j})  -label "$j " \
		-step 1 -selectmode immediate -options { entry.width 3 }
	pack $w0.evg${i}_${j} -in $fmri(modelf$i) -side top -anchor w -padx 5 -pady 2
    }

}

#}}}
		} else {
		    $w0.evsnb pageconfigure ev$i -state normal
		}

		feat5:setup_model_update_ev_i $w $w0 $i 0 0

	    } else {
		#{{{ disable nb page

		if { [ info exists fmri(modelf$i) ] && [ winfo exists $fmri(modelf$i) ] } {
		    $w0.evsnb pageconfigure ev$i -state disabled
		} else {
		    set i -10
		}

#}}}
	    }
	}

#}}}
	} else {
	    #{{{ EV grid for level>1

set w1 $fmri(evsf).grid
if { [ winfo exists $w1 ] } {
    destroy $w1
}
frame $w1
pack $w1 -in $fmri(evsf) -side top -anchor w -padx 5 -pady 5

button $w1.pastebutton -command "feat5:multiple_paste \"Higher-level model EVs\" $fmri(evs_orig) $fmri(multiple) fmri evg" -text "Paste"
grid $w1.pastebutton -in $w1 -row 0 -column 0

label $w1.grouplabel -text "     Group     "
grid $w1.grouplabel -in $w1 -row 0 -column 1

for { set j 1 } { $j <= $fmri(evs_orig) } { incr j 1 } {
    label $w1.evlabel${j} -text "EV$j"
    grid $w1.evlabel${j} -in $w1 -row 0 -column [ expr $j + 1 ]
    entry $w1.evtitle${j} -textvariable fmri(evtitle$j) -width 7
    grid $w1.evtitle${j} -in $w1 -row 1 -column [ expr $j + 1 ] -padx 2 -pady 2
    $w.bhelp bind $w1.evtitle${j} -msg "If you wish, enter a title for EV $j here."
}

for { set i 1 } { $i <= $fmri(multiple) } { incr i 1 } {

    label $w1.label$i -text "Input $i "
    grid $w1.label$i -in $w1 -row [ expr $i + 1 ] -column 0
    $w.bhelp bind $w1.label$i -msg "Input (lower-level FEAT directory) ${i}."

    if { ! [ info exists fmri(groupmem.$i) ] } {
	set fmri(groupmem.$i) 1
    }
    tixControl $w1.groupmem$i -variable fmri(groupmem.$i) -step 1 -min 1 -selectmode immediate -options { entry.width 3 }
    grid $w1.groupmem$i -in $w1 -row [ expr $i + 1 ] -column 1
    $w.bhelp bind $w1.groupmem$i -msg "Which group of subjects (or sessions, etc.) is this input a part of?

If you setup different groups for different variances, you will get
fewer data-points to estimate each variance (than if only one variance
was estimated). Therefore, you only want to use this option if you do
believe that the groups possibly do have different variances.

If you setup different groups for different variances, it is necessary
that, for each EV, only one of the sub-groups has non-zero
values. Thus, for example, in the case of an unpaired t-test:

GP EV1 EV2
1   1    1 
1   1    1 
1   1    1 
1   1    1 
2   1   -1
2   1   -1
2   1   -1

is wrong with respect to this issue, and the following is correct:

GP EV1 EV2
1   1    0 
1   1    0 
1   1    0 
1   1    0 
2   0    1
2   0    1
2   0    1"

    for { set j 1 } { $j <= $fmri(evs_orig) } { incr j 1 } {
	tixControl $w1.evg${i}_${j} -variable fmri(evg${i}.${j}) -step 1 -selectmode immediate -options { entry.width 3 }
	grid $w1.evg${i}_${j} -in $w1 -row [ expr $i + 1 ] -column [ expr $j + 1 ]
	$w.bhelp bind $w1.evg${i}_${j} -msg "Design matrix value for Input (lower-level FEAT directory) $i and EV ${j}."
    }
}

if { [ winfo exists $fmri(evsf).orthbutton ] } {
    destroy $fmri(evsf).orthbutton
}
if { $fmri(evs_orig) > 1 } {
    if { ! [ info exists fmri(level2orth) ] || ! $fmri(level2orth) } {
	button $fmri(evsf).orthbutton -command "feat5:setup_level2orth $w $w1" -text "Setup orthogonalisations"
	pack $fmri(evsf).orthbutton -in $fmri(evsf) -side top -anchor w -padx 5 -pady 5
	set fmri(level2orth) 0
    } else {
	feat5:setup_level2orth $w $w1
    }
}

#}}}
	}
    }

    feat5:setup_model_update_contrasts $w 0
}    

#}}}
#{{{ feat5:setup_model_update_contrasts_real_per_orig

proc feat5:setup_model_update_contrasts_real_per_orig { w } {
    global fmri

    # first - are there any basis functions OR sinusoidal harmonics? if not return 0
    set is_simple 1
    for { set i 1 } { $i <= $fmri(evs_orig) } { incr i 1 } {
	if { $fmri(convolve$i) > 3 || $fmri(shape$i) == 1 } {
	    set is_simple 0
	}
    }
    if { $is_simple } {
	return 0
    }

    # if there are make sure the number of BF or harmonics is the same for every original EV
    for { set i 2 } { $i <= $fmri(evs_orig) } { incr i 1 } {
	if { $fmri(shape$i) != $fmri(shape1) || $fmri(evs_real.$i) != $fmri(evs_real.1) } {
	    return -1
	}
    }

    # so we're ok - return the number of BF or harmonics
    return $fmri(evs_real.1)	
}

#}}}
#{{{ feat5:setup_model_update_contrasts_mode

proc feat5:setup_model_update_contrasts_mode { w update_gui } {
    global fmri

    # if called from feat5:write, will have been called with mode=C
    if { ! $update_gui } {
	set fmri(con_mode_old) c
    }

    if { $fmri(con_mode) != $fmri(con_mode_old) } {

	if { $update_gui } {
	    $fmri(contrastsf).num.con    configure -variable fmri(ncon_$fmri(con_mode))
	    $fmri(contrastsf).num.ftests configure -variable fmri(nftests_$fmri(con_mode))
	}

	if { $fmri(con_mode) == "real" || ! $update_gui } {

	    set real_per_orig [ feat5:setup_model_update_contrasts_real_per_orig $w ]
	    if { $real_per_orig == -1 } {
		MxPause "In order to setup contrasts in \"Original EVs\" mode whilst using basis functions or sinusoidal harmonics, all the original EVs must be the same basic shape and generate the same number of real EVs. Please change the EV setup."
		return -1
	    }

	    if { $real_per_orig > 0 } {
		#{{{ do the case of basis functions or sinusoidal harmonics

set fmri(ncon_real)    [ expr $real_per_orig * $fmri(ncon_orig) ]
set fmri(nftests_real) [ expr $fmri(nftests_orig) + $fmri(ncon_orig) ]

# zero all F-tests
for { set Con 1 } { $Con <= $fmri(ncon_real) } { incr Con 1 } {
    for { set F 1 } { $F <= $fmri(nftests_real) } { incr F 1 } {
	set fmri(ftest_real${F}.$Con) 0
    }
}

# set explicitly asked for F-tests
for { set F 1 } { $F <= $fmri(nftests_orig) } { incr F 1 } {
    for { set Con 1 } { $Con <= $fmri(ncon_orig) } { incr Con 1 } {
	if { $fmri(ftest_orig${F}.$Con) == 1 } {
	    for { set con_real_inc 1 } { $con_real_inc <= $real_per_orig } { incr con_real_inc } {
		set fmri(ftest_real${F}.[ expr ( ( $Con - 1 ) * $real_per_orig ) + $con_real_inc ]) 1
	    }
	}
    }
}

# set the rest
set con_real 0
for { set Con 1 } { $Con <= $fmri(ncon_orig) } { incr Con 1 } {

    set ev_real 0
    for { set ev_orig 1 } { $ev_orig <= $fmri(evs_orig) } { incr ev_orig 1 } {
	for { set con_real_inc 1 } { $con_real_inc <= $real_per_orig } { incr con_real_inc } {
	    set fmri(conpic_real.[ expr $con_real + $con_real_inc ]) $fmri(conpic_orig.$Con)
	    set fmri(conname_real.[ expr $con_real + $con_real_inc ]) "$fmri(conname_orig.$Con) ($con_real_inc)"
	    for { set ev_real_inc 1 } { $ev_real_inc <= $real_per_orig } { incr ev_real_inc } {
		if { $con_real_inc == $ev_real_inc } {
		    set fmri(con_real[ expr $con_real + $con_real_inc ].[ expr $ev_real + $ev_real_inc ]) $fmri(con_orig${Con}.$ev_orig)
		} else {
		    set fmri(con_real[ expr $con_real + $con_real_inc ].[ expr $ev_real + $ev_real_inc ]) 0
		}
	    }
	}
	incr ev_real $real_per_orig
    }

    for { set con_real_inc 1 } { $con_real_inc <= $real_per_orig } { incr con_real_inc } {
	set fmri(ftest_real[ expr $Con + $fmri(nftests_orig) ].[ expr $con_real + $con_real_inc ]) 1
    }

    incr con_real $real_per_orig
}

#}}}
	    } else {
		#{{{ do temporal derivates etc.

for { set Con 1 } { $Con <= $fmri(ncon_orig) } { incr Con 1 } {
    set ev_real 1
    for { set ev_orig 1 } { $ev_orig <= $fmri(evs_orig) } { incr ev_orig 1 } {
	set fmri(conpic_real.$Con) $fmri(conpic_orig.$Con)
	set fmri(conname_real.$Con) $fmri(conname_orig.$Con)
	set fmri(con_real${Con}.$ev_real) $fmri(con_orig${Con}.$ev_orig)
	incr ev_real 1
	if { $fmri(deriv_yn$ev_orig) } {
	    set fmri(con_real${Con}.$ev_real) 0
	    incr ev_real 1
	}
    }    

    for { set F 1 } { $F <= $fmri(nftests_orig) } { incr F 1 } {
	set fmri(ftest_real${F}.${Con}) $fmri(ftest_orig${F}.${Con})
    }
}

set fmri(ncon_real) $fmri(ncon_orig)
set fmri(nftests_real) $fmri(nftests_orig)

#}}}
	    }

	}

	set fmri(con_mode_old) $fmri(con_mode)
    }

    return 1
}

#}}}
#{{{ feat5:setup_model_update_contrasts

proc feat5:setup_model_update_contrasts { w dummy } {
    global fmri

    #{{{ setup evs_real etc/

    set fmri(evs_real) 0

    for { set j 1 } { $j <= $fmri(evs_orig) } { incr j 1 } {
	set fmri(evs_real.$j) 1

	incr fmri(evs_real.$j) $fmri(deriv_yn$j)

	if { $fmri(convolve$j) > 3 } {
	    incr fmri(evs_real.$j) [ expr $fmri(basisfnum$j) - 1 ]
	}

	if { $fmri(shape$j) == 1 } {
	    incr fmri(evs_real.$j) $fmri(nharmonics$j)
	}

	incr fmri(evs_real) $fmri(evs_real.$j)
    }

#}}}

    if { $fmri(level) > 1 } {
	set fmri(con_mode) real
    }

    for { set i 1 } { $i <= $fmri(ncon_$fmri(con_mode)) } { incr i 1 } {
	if { ! [ info exists fmri(conpic_$fmri(con_mode).$i) ] } {
	    set fmri(conpic_$fmri(con_mode).$i) 1
	    set fmri(conname_$fmri(con_mode).$i) ""
	}
    }

    #{{{ destroy and recreate grid

set w0 $fmri(contrastsf).congrid

if { [ winfo exists $w0 ] } {
    destroy $w0
}

frame $w0

pack $w0 -in $fmri(contrastsf) -side top -anchor w -padx 5 -pady 5

#}}}
    #{{{ first 3 columns

button $w0.pastebutton -command "feat5:multiple_paste \"Contrasts\" $fmri(evs_$fmri(con_mode)) $fmri(ncon_$fmri(con_mode)) fmri con_$fmri(con_mode)" -text "Paste"
grid $w0.pastebutton -in $w0 -row 0 -column 0

for { set i 1 } { $i <= $fmri(ncon_$fmri(con_mode)) } { incr i 1 } {

    if { $fmri(con_mode) == "orig" } {
	label $w0.label$i -text "OC$i "
    } else {
	label $w0.label$i -text "C$i "
    }

    grid $w0.label$i -in $w0 -row $i -column 0
    $w.bhelp bind $w0.label$i -msg "Contrast vector number $i - this will result in Z statistic image number $i"

    checkbutton $w0.conpic$i -variable fmri(conpic_$fmri(con_mode).$i)
    grid $w0.conpic$i -in $w0 -row $i -column 1
    $w.bhelp bind $w0.conpic$i -msg "Include contrast $i in web page report? (Turn off if
this contrast is only to be used within F-tests.)"

    entry $w0.conname$i -textvariable fmri(conname_$fmri(con_mode).$i) -width 10
    grid $w0.conname$i -in $w0 -row $i -column 2
}

#}}}
    #{{{ top row and main grid

label $w0.evlabel0 -text "Title"
grid $w0.evlabel0 -in $w0 -row 0 -column 2

set ev_count 1
set ev_counti 1

for { set j 1 } { $j <= $fmri(evs_$fmri(con_mode)) } { incr j 1 } {

    if { $fmri(con_mode) == "real" || $ev_counti == 1 } {

	if { $ev_counti == 1 } {
	    label $w0.evlabel${j} -text "EV$ev_count"
	} else {
	    label $w0.evlabel${j} -text ""
	}
	grid $w0.evlabel${j} -in $w0 -row 0 -column [ expr $j + 2 ]

	for { set i 1 } { $i <= $fmri(ncon_$fmri(con_mode)) } { incr i 1 } {
	    if { ! [ info exists fmri(con_$fmri(con_mode)${i}.${j}) ] } {
		set fmri(con_$fmri(con_mode)${i}.${j}) 0
		if { $i==1 && $j==1 } {
		    set fmri(con_$fmri(con_mode)${i}.${j}) 1
		}
	    }
	    tixControl $w0.con${i}_${j} -variable fmri(con_$fmri(con_mode)${i}.${j}) -step 1 -selectmode immediate -options { entry.width 3 }
	    grid $w0.con${i}_${j} -in $w0 -row $i -column [ expr $j + 2 ]
	    $w.bhelp bind $w0.con${i}_${j} -msg "The weight given to EV$j within contrast vector $i."
	}


    }

    if { $fmri(con_mode) == "real" } {
	incr ev_counti 1
	if { $ev_counti > $fmri(evs_real.$ev_count) } {
	    set ev_counti 1
	    incr ev_count 1
	}
    } else {
	incr ev_count 1
    }
}

#}}}
    #{{{ ftests

label $w0.blank -text "       "
grid $w0.blank -in $w0 -row 0 -column [ expr $fmri(evs_real) + 3 ]

for { set i 1 } { $i <= $fmri(nftests_$fmri(con_mode)) } { incr i 1 } {

    label $w0.flabel$i -text "F$i "
    grid $w0.flabel$i -in $w0 -row 0 -column [ expr $i + $fmri(evs_real) + 3 ]
    $w.bhelp bind $w0.flabel$i -msg "F-test number $i - this will result in F statistic image number $i"
	
    for { set j 1 } { $j <= $fmri(ncon_$fmri(con_mode)) } { incr j 1 } {
	
	if { ! [ info exists fmri(ftest_$fmri(con_mode)${i}.${j}) ] } {
	    set fmri(ftest_$fmri(con_mode)${i}.${j}) 0
	}

	checkbutton $w0.ftest${i}_${j} -variable fmri(ftest_$fmri(con_mode)${i}.$j)
	grid $w0.ftest${i}_${j} -in $w0 -row $j -column [ expr $i + $fmri(evs_real) + 3 ]
	$w.bhelp bind $w0.ftest${i}_${j} -msg "Include contrast vector $j in F-test $i?"
    }
}

#}}}
}

#}}}
#{{{ feat5:setup_model_preview

proc feat5:setup_model_preview { w } {
    global fmri

    set fmri(filmsetup) 1

    set problem [ feat5:write $w 1 0 0 $fmri(feat_filename) ]
    if { $problem } {
	return 1
    }

    set count 0
    set w0 ".dialog[incr count]"
    while { [ winfo exists $w0 ] } {
        set w0 ".dialog[incr count]"
    }

    toplevel $w0 -visual truecolor
    wm title $w0 "Model"
    wm iconname $w0 "Model"

    frame $w0.f
    pack $w0.f -expand yes -fill both -in $w0 -side top
    canvas $w0.f.viewport -yscrollcommand "$w0.f.ysbar set" -xscrollcommand "$w0.xsbar set"
    scrollbar $w0.xsbar -command "$w0.f.viewport xview" -orient horizontal
    scrollbar $w0.f.ysbar -command "$w0.f.viewport yview" -orient vertical
    frame $w0.f.viewport.f
    $w0.f.viewport create window 0 0 -anchor nw -window $w0.f.viewport.f
    bind $w0.f.viewport.f <Configure> "feat5:scrollform_resize $w0 $w0.f.viewport"
    pack $w0.f.viewport -side left -fill both -expand true -in $w0.f

    set graphpic [ image create photo -file [ file rootname $fmri(feat_filename) ].ppm ]
    button $w0.f.viewport.f.btn -command "destroy $w0" -image $graphpic -borderwidth 0
    pack $w0.f.viewport.f.btn -in $w0.f.viewport.f

    $w.bhelp bind $w0.f.viewport.f.btn -msg $fmri(design_help)

    return 0
}

#}}}
#{{{ feat5:setup_model_acpreview

proc feat5:setup_model_acpreview { w } {
    global fmri

    set problem [ feat5:write $w 1 0 0 $fmri(feat_filename) ]
    if { $problem } {
	return 1
    }

    set count 0
    set w1 ".dialog[incr count]"
    while { [ winfo exists $w1 ] } {
        set w1 ".dialog[incr count]"
    }

    toplevel $w1 -visual truecolor
    wm title $w1 "Design efficiency"
    wm iconname $w1 "Efficiency"

    set graphpic [ image create photo -file [ file rootname $fmri(feat_filename) ]_cov.ppm ]

    button $w1.btn -command "destroy $w1" -image $graphpic -borderwidth 0
    pack $w1.btn -in $w1

    $w.bhelp bind $w1.btn -msg $fmri(cov_help)

    return 0
}

#}}}
#{{{ feat5:setup_model_destroy

proc feat5:setup_model_destroy { w w0 } {
    global fmri

    if { ! [ feat5:setup_model_preview $w ] } {
	destroy $w0
    }
}

#}}}
#{{{ feat5:setup_conmask

proc feat5:setup_conmask { w } {
    global fmri FSLDIR

    #{{{ setup window

    if { [ info exists fmri(c_model) ] && [ winfo exists $fmri(c_model) ] } {
	return 0
    }

    set count 0
    set w0 ".cdialog[incr count]"
    while { [ winfo exists $w0 ] } {
        set w0 ".cdialog[incr count]"
    }

    set fmri(c_model) $w0

    toplevel $w0

    wm title $w0 "Setup Contrast Masking"
    wm iconname $w0 "Contrast Masking"
    wm iconbitmap $w0 @${FSLDIR}/tcl/fmrib.xbm

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
    #{{{ setup grid

set total [ expr $fmri(ncon_real) + $fmri(nftests_real) ]

for { set C 1 } { $C <= $total } { incr C 1 } {
    if { $C <= $fmri(ncon_real) } {
	label $v.tl$C -text "C$C"
    } else {
	label $v.tl$C -text "F[ expr $C - $fmri(ncon_real)]"
    }
    $w.bhelp bind $v.tl$C -msg $fmri(conmask_help) 
    grid $v.tl$C -in $v -column $C -row 0
}

for { set c 1 } { $c <= $total } { incr c 1 } {

    if { $c <= $fmri(ncon_real) } {
	label $v.l$c -text "Mask real Contrast $c with:    "
    } else {
	label $v.l$c -text "Mask real F-test [ expr $c - $fmri(ncon_real)] with:    "
    }
    $w.bhelp bind $v.l$c -msg $fmri(conmask_help) 
    grid $v.l$c -in $v -column 0 -row $c

    for { set C 1 } { $C <= $total } { incr C 1 } {

	if { $C != $c } {
	    checkbutton $v.cb${c}_$C -variable fmri(conmask${c}_${C})
	    $w.bhelp bind $v.cb${c}_$C -msg $fmri(conmask_help) 
	    grid $v.cb${c}_$C -in $v -column $C -row $c
	}
    }
}

#}}}
    #{{{ setup buttons

checkbutton $w0.zeros -variable fmri(conmask_zerothresh_yn) -text "Mask using (Z>0) instead of (Z stats pass thresholding)" 
$w.bhelp bind $w0.zeros -msg $fmri(conmask_help) 

button $w0.ok -command "destroy $w0" -text "OK"

pack $w0.ok $w0.zeros -in $w0 -side bottom -padx 3 -pady 5

#}}}
}

#}}}
#{{{ feat5:checkbfcustom

proc feat5:checkbfcustom { w i dummy } {
    global fmri

    if { [ string compare [ file extension $fmri(bfcustom$i) ] .flobs ] == 0 } {
	set fmri(bfcustom$i) $fmri(bfcustom$i)/hrfbasisfns.txt
    }

    if { ! [ file exists $fmri(bfcustom$i) ] } {
	MxPause "Custom HRF convolution file is invalid!"
	return 1
    }

    catch { exec sh -c "wc -l $fmri(bfcustom$i) | awk '{ print \$1 }'" } line_count
    catch { exec sh -c "wc -w $fmri(bfcustom$i) | awk '{ print \$1 }'" } word_count

    set fmri(basisfnum$i) [ expr int ( $word_count / $line_count ) ]
}

#}}}
#{{{ feat5:setup_level2orth

proc feat5:setup_level2orth { w w1 } {
    global fmri

    set fmri(level2orth) 1

    destroy $fmri(evsf).orthbutton $fmri(evsf).orthlabel 

    label $fmri(evsf).orthlabel -text "Orthogonalisations  "
    grid $fmri(evsf).orthlabel -in $w1 -row [ expr 2 + $fmri(multiple) ] -column 0 -columnspan 2

    for { set i 1 } { $i <= $fmri(evs_orig) } { incr i 1 } {
	set fmri(ortho${i}.0) 1
	set fmri(ortho${i}.${i}) 0
	for { set j 1 } { $j <= $fmri(evs_orig) } { incr j 1 } {
	    if { $j != $i } {
		checkbutton $w1.doorth${i}_${j} -variable fmri(ortho${i}.$j)
		grid $w1.doorth${i}_${j} -in $w1 -row [ expr $j + $fmri(multiple) + 1 ] -column [ expr $i + 1 ]
		$w.bhelp bind $w1.doorth${i}_${j} -msg "Orthogonalise EV$i wrt EV${j}?"
	    }
	}
    }
}

#}}}

#}}}
#{{{ feat5:estnoise

proc feat5:estnoise { } {

    global fmri feat_files FSLDIR

    set smooth [ expr $fmri(smooth) / 2.355 ]

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

    if { [ info exists feat_files(1) ] && [ imtest $feat_files(1) ] } {
	set noiseestout [ exec sh -c "${FSLDIR}/bin/estnoise $feat_files(1) $smooth $hp_sigma_vol $lp_sigma_vol 2> /dev/null" ]
	set fmri(noise)   [ lindex $noiseestout 0 ]
	set fmri(noisear) [ lindex $noiseestout 1 ]
    }
}

#}}}
