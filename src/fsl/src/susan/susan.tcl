#

#{{{ setup

#   susan.tcl - GUI for susan smoothing
#
#   Stephen Smith, FMRIB Image Analysis Group
#
#   Copyright (C) 1999-2004 University of Oxford
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

#{{{ proc susan

proc susan { w } {

    global susanvars usanentries FSLDIR PWD argc argv TN HOME
 
    #{{{ setup main window

toplevel $w

wm title $w "SUSAN   Structure-Preserving Noise Reduction"
wm iconname $w "SUSAN"
wm iconbitmap $w @$FSLDIR/tcl/fmrib.xbm

frame $w.f

#}}}
    #{{{ setup variable defaults

set susanvars($w,dt)  3
set susanvars($w,um)  1
set susanvars($w,dim) 3D

#}}}
    #{{{ input and output images

    if { $argc > 0 && [ string length [ lindex $argv 0 ] ] > 0 } {
	set inputname [ remove_ext [ imglob -oneperimage [ lindex $argv 0 ] ] ]
        if { [ imtest $inputname ] } {
	    if { [ string first / $inputname ] == 0 || [ string first ~ $inputname ] == 0 } {
		set susanvars($w,input) $inputname
	    } else {
		set susanvars($w,input) ${PWD}/$inputname
	    }
	    set susanvars($w,output) $susanvars($w,input)_susan
	    susan:setrange $w dummy
	}
    }

    FSLFileEntry $w.f.input \
	    -variable susanvars($w,input) \
	    -pattern "IMAGE" \
	    -directory $PWD \
	    -label "Input image   " \
	    -title "Select the input image" \
	    -width 55 \
	    -command "susan:setrange $w" \
	    -filterhist VARS(history)

    FSLFileEntry $w.f.output \
	    -variable susanvars($w,output) \
	    -pattern "IMAGE" \
	    -directory $PWD \
	    -label "Output image" \
	    -title "Select the output image" \
	    -width 55 \
	    -filterhist VARS(history)

    pack $w.f.input $w.f.output -in $w.f -side top -padx 5 -pady 5 -anchor w


#}}}
    #{{{ top row

    frame $w.f.top

    tixOptionMenu $w.f.op -label "Dimensionality" 
    $w.f.op add command 2D -label 2D
    $w.f.op add command 3D -label 3D
    $w.f.op config -variable susanvars($w,dim) 

    tixControl $w.f.bt -label "Brightness threshold" -variable susanvars($w,bt) -step 1 -min 0 -selectmode immediate -options {
	entry.width 8
    }

    tixControl $w.f.dt -label "Mask SD" -variable susanvars($w,dt) -step 1 -min 0 -selectmode immediate -options {
	entry.width 5
    }

    pack $w.f.top -in $w.f -side top -padx 3 -pady 3 -expand yes -anchor w

    pack $w.f.op $w.f.bt $w.f.dt -in $w.f.top -side left -padx 5 -pady 3 -expand yes -anchor w

#}}}
    #{{{ help text

label $w.f.label -wraplength 550 -text "SUSAN filtering works on 16 bit signed integer data.\nMask SD (approx HWHM) is in current length units (e.g. mm), not voxels.\nSet Mask SD to 0 for a (much faster) flat 3x3 voxels or 3x3x3 voxels mask." -fg "#303030"

pack $w.f.label -in $w.f -side top -anchor s

#}}}
    #{{{ advanced options

collapsible frame $w.f.opts -title "Advanced options"

#{{{ use median?

frame $w.f.um

label $w.f.um.label -text "Use median when no neighbourhood is found "

checkbutton $w.f.um.but -variable susanvars($w,um)

pack $w.f.um.label $w.f.um.but -in $w.f.um -side left

pack $w.f.um  -in $w.f.opts.b -anchor w -pady 5

#}}}
#{{{ more USANs

set susanvars($w,maxusans) 2
set susanvars($w,usans) 0

tixControl $w.f.usans -label "Separate images to find USAN from " \
	-variable susanvars($w,usans) -step 1 -min 0 -max $susanvars($w,maxusans) -selectmode immediate -command "susan:updateusan $w"

pack $w.f.usans -in $w.f.opts.b -anchor w -pady 5

set i 1
while { $i <= $susanvars($w,maxusans) } {
    frame $w.f.usanentries($i)

	FSLFileEntry $w.f.ue$i \
		-variable usanentries($w,$i) \
		-pattern "IMAGE" \
		-directory $PWD \
		-label "USAN image $i" \
		-title "Select the USAN image" \
		-width 40 \
		-filterhist VARS(history)
	
    set susanvars($w,ubt,$i) 1
    tixControl $w.f.ubt$i -label "Brightness threshold" -variable susanvars($w,ubt,$i) -step 1 -min 0 -selectmode immediate -options {
	entry.width 5
    }

    pack $w.f.ue$i $w.f.ubt$i -in $w.f.usanentries($i) -padx 3 -pady 3 -side left
    incr i 1
}

#}}}

pack $w.f.opts -in $w.f -side bottom -anchor w -pady 5

#}}}
    #{{{ button frame

    frame $w.btns
    frame $w.btns.b -relief raised -borderwidth 1
    
    button $w.ok -text "OK" -width 5 -command "susan:apply $w destroy"
    bind $w.ok <Return> {
        [winfo toplevel %W].ok invoke
    }
    
    button $w.apply -command "susan:apply $w keep" -text "Apply" -width 5
    bind $w.apply <Return> {
	[winfo toplevel %W].apply invoke
    }
	    
    button $w.cancel    -command "susan:destroy $w" \
	    -text "Cancel" -width 5
    bind $w.cancel <Return> {
	[winfo toplevel %W].cancel invoke
    }

    button $w.help -command "FmribWebHelp file: $FSLDIR/doc/susan/index.html" \
            -text "Help" -width 5
    bind $w.help <Return> {
        [winfo toplevel %W].help invoke
    }

    pack $w.btns.b -side bottom -fill x -padx 3 -pady 3
    pack $w.ok $w.apply $w.cancel $w.help -in $w.btns.b \
	    -side left -expand yes -padx 3 -pady 3 -fill y

#}}}
		
    pack $w.f $w.btns -expand yes -fill both
}

#}}}
#{{{ proc susan:setrange

proc susan:setrange { w dummy } {

    global susanvars FSLDIR

    set susanvars($w,input)  [ remove_ext $susanvars($w,input) ]
    set susanvars($w,output) [ remove_ext $susanvars($w,input) ]_susan
    
    if { ! [ catch { exec sh -c "${FSLDIR}/bin/avwstats $susanvars($w,input) -r" } minmax ] } {
	set min [ lindex $minmax 0 ]
	set max [ lindex $minmax 1 ]
	set susanvars($w,bt) [ expr ( $max - $min ) / 10.0 ]
    }
}

#}}}
#{{{ proc susan:apply

proc susan:apply { w dialog } {

    global susanvars usanentries TN

    foreach v { bt dt } {
	$w.f.$v update
    }

    #{{{ process USAN entries

set usanlist ""

set susanvars($w,name) grot

set i 1
while { $i <= $susanvars($w,usans) } {
    lappend usanlist [ remove_ext $usanentries($w,$i) ]
    lappend usanlist $susanvars($w,ubt,$i)
    incr i 1
}

#}}}
    
    susan_proc $susanvars($w,input) $susanvars($w,name) $susanvars($w,bt) $susanvars($w,output) $susanvars($w,dt) $susanvars($w,dim) $susanvars($w,um) $susanvars($w,usans) $usanlist
    
    update idletasks
    
    if {$dialog == "destroy"} {
	susan:destroy $w
    }
}

#}}}
#{{{ proc susan:destroy

proc susan:destroy { w } {
    destroy $w
}

#}}}
#{{{ updates

proc susan:updateusan { w val } {
    global susanvars

    if { $susanvars($w,usans) == 0 } {
	pack $w.f.bt -in $w.f.top -before $w.f.dt -side left -padx 5 -pady 3 -expand yes -anchor w
    } else {
	pack forget $w.f.bt
    }

    set i 1
    while { $i <= $susanvars($w,maxusans) } {
	pack forget $w.f.usanentries($i)
	incr i 1
    }

    set i 1
    while { $i <= $susanvars($w,usans) } {
	pack $w.f.usanentries($i) -in $w.f.opts.b -anchor w -pady 5
	incr i 1
    }
}

#}}}

#{{{ tail end

wm withdraw .
susan .rename
tkwait window .rename

#}}}
