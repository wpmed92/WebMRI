#   FSL interface for FDT (BEDPOST and ProbTrack)
#
#   Timothy Behrens, Heidi Johansen-Berg and Dave Flitney, FMRIB Image Analysis Group
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
set TCLPATH [file dirname [ info script ] ]
regsub tcl $TCLPATH bin BINPATH
regsub tcl $TCLPATH doc/fdt HTMLPATH

set VERSION "1.0a"

tixWidgetClass registrationImageSelect {
    -superclass tixPrimitive
    -classname RegistrationImageSelect
    
    -flag {
	-variable -filename -labelwidth -label -search -dof -pattern -filterhist -costfn
    }
    
    -configspec {
	{-variable variable Variable N}
	{-labelwidth labelwidth Labelwidth {}}
	{-label label Label ""}
	{-directory directory Directory {}}
	{-filterhist filterhist Filterhist {}}
	{-filename filename Filename {}}
	{-search search Search N}
	{-dof dof DOF N}
	{-costfn costfn Costfn N}
	{-pattern pattern Pattern "*.{hdr,hdr.gz,nii,nii.gz}"}
    }
}

proc registrationImageSelect:ConstructWidget { w } {

    upvar \#0 $w data

    tixChainMethod $w ConstructWidget

    set data(w:fsd) [ tixFileSelectDialog $w.fsd ]

    frame $w.frame
    
    set data(w:chkbutton) [ checkbutton $w.frame.visible_yn -variable $data(-variable) \
				-command "registrationImageSelect:toggle_visible $w" ]
    label $w.frame.label -text $data(-label)

    set data(w:options) [ tixLabelFrame $w.frame.options -label $data(-label) ]

    set data(w:file) [ FSLFileEntry $w.frame.options.f \
			   -variable $data(-filename) \
			   -directory $data(-directory) \
			   -label "" \
			   -title "Select file" \
			   -width 40 \
			   -pattern $data(-pattern) \
			   -filterhist VARS(history) ]

    tixOptionMenu $w.frame.options.search -variable $data(-search)
    $w.frame.options.search add command 0   -label "No search"
    $w.frame.options.search add command 90  -label "Normal search"
    $w.frame.options.search add command 180 -label "Full search"

    tixOptionMenu $w.frame.options.dof -variable $data(-dof)
    $w.frame.options.dof add command 6  -label "6 DOF"
    $w.frame.options.dof add command 7  -label "7 DOF"
    $w.frame.options.dof add command 9  -label "9 DOF"
    $w.frame.options.dof add command 12 -label "12 DOF"

    tixOptionMenu $w.frame.options.costfn -variable $data(-costfn)
    $w.frame.options.costfn add command corratio -label "Correlation ratio"
    $w.frame.options.costfn add command mutualinfo -label "Mutual information"
    
    pack $w.frame.options.f -in [$w.frame.options subwidget frame] -side top -expand yes -fill x
    pack $w.frame.options.costfn $w.frame.options.dof $w.frame.options.search -in [$w.frame.options subwidget frame] -side right
    pack $w.frame.visible_yn $w.frame.label -in $w.frame -side left -expand yes -padx 3
    pack $w.frame -in $w
    
    registrationImageSelect:toggle_visible $w
}

proc registrationImageSelect:toggle_visible { w } {

    upvar \#0 $w data
    upvar \#0 $data(-variable) visible_yn

    if { $visible_yn } {
	pack forget $w.frame.label
	pack [$w subwidget options] -side left -after [$w subwidget chkbutton] -padx 3 -pady 3
    } else {
	pack $w.frame.label -side left -after [$w subwidget chkbutton] -padx 3 -pady 3
	pack forget [$w subwidget options]
    }
}

tixWidgetClass multiFileSelect {
    -superclass tixPrimitive
    -classname MultiFileSelect
    -method {
	save load
    }
    -flag {
	-variable -labelwidth -label -directory -pattern -title -filterhist
    }

    -configspec {
	{-variable variable Variable N}
	{-labelwidth labelwidth Labelwidth {}}
	{-label label Label "File"}
	{-directory directory Directory {}}
	{-filterhist filterhist Filterhist {}}
	{-pattern pattern Pattern "*.{hdr,hdr.gz,nii,nii.gz}"}
	{-title title Title "Select a file"}
    }
}

proc multiFileSelect:InitWidgetRec { w } {

    upvar #0 $w data

    tixChainMethod $w InitWidgetRec

    set data(filename) filename$w
}

proc multiFileSelect:ConstructWidget { w } {

    upvar #0 $w data

    tixChainMethod $w ConstructWidget

    set data(w:fsd) [ tixFileSelectDialog $w.fsd ]

    frame $w.frame -borderwidth 2 -relief sunken

    frame $w.frame.l
    set data(w:files) [ listbox $w.frame.l.files -height 4 -yscrollcommand "$w.frame.l.scroller set" ]

    scrollbar $w.frame.l.scroller -command "$w.frame.l.files yview"
    pack $w.frame.l.files -in $w.frame.l -side left -expand yes -fill x
    pack $w.frame.l.scroller -in $w.frame.l -side right -fill y
 
    set data(w:file) [FSLFileDirSelectDialog $w.browser \
	    -title $data(-title) \
	    -command "multiFileSelect:Invoke $w" \
	    -filterhist $data(-filterhist) \
	    -directory $data(-directory) \
	    -pattern $data(-pattern)
    ]

    frame $w.frame.b

    button $w.frame.b.add -text "Add..."  -command "$w.browser popup"
    button $w.frame.b.del -text "Remove"  -command "multiFileSelect:remove $w"
    button $w.frame.b.imp -text "Load list..." -command "multiFileSelect:import $w"
    button $w.frame.b.exp -text "Save list..." -command "multiFileSelect:export $w"

    pack $w.frame.b.add $w.frame.b.del $w.frame.b.imp $w.frame.b.exp -in $w.frame.b -side left -expand yes -fill x
 
    pack $w.frame.l $w.frame.b -in $w.frame -side top -padx 2 -pady 2 -expand yes -fill x

    pack $w.frame -in $w -expand yes -fill x 
}

proc multiFileSelect:Invoke { w filename } {
    upvar #0 $w data
    upvar $data(-variable) Variable
    
    set filename [ fix_cygwin_filename $filename ]
    set Variable $filename

    [$w subwidget files] insert end $filename
}

proc multiFileSelect:remove { w } {
    upvar #0 $w data
    set count 0
    foreach file [ [$w subwidget files] get 0 end ] {
	if { [ [$w subwidget files] selection includes $count ] == 1 } {
	    [$w subwidget files] delete $count
	    incr count -1
	}
	incr count
    } 
}

proc multiFileSelect:load { w filename } {
    upvar #0 $w data

    if { ![ file exists $filename ] } {
	MxPause "Warning: Bad or missing file!"
	return
    }
    set fd [ open $filename ]
    [$w subwidget files] delete 0 end
    while { [ gets $fd file ] >= 0 } {
	[$w subwidget files] insert end $file
    }
    close $fd
}

proc multiFileSelect:save { w filename } {
    upvar #0 $w data

    set fd [ open $filename w ]
    foreach file [ [$w subwidget files] get 0 end ] {
	puts $fd $file
    }
    close $fd
}

proc multiFileSelect:import { w } {   
    upvar #0 $w data

    set fsd $data(w:fsd)
    $fsd configure -title "Import list from..." -command "multiFileSelect:load $w"
    $fsd popup
}

proc multiFileSelect:export { w } {
    upvar #0 $w data

    set fsd $data(w:fsd)
    $fsd configure -title "Export list to..." -command "multiFileSelect:save $w"
    $fsd popup
}

tixWidgetClass coordinateEdit {
    -superclass tixPrimitive
    -classname CoordinateEdit
    -method {
	update
    }
    -flag {
	-variablex -variabley -variablez -variablemm -labelwidth -label
    }

    -configspec {
	{-variablex variablex VariableX N}
	{-variabley variabley VariableY N}
	{-variablez variablez VariableZ N}
	{-variableunits variableunits VariableUnits N}
	{-label label Label "File"}
	{-labelwidth labelwidth Labelwidth {}}
    }
}

proc coordinateEdit:ConstructWidget { w } {

    upvar #0 $w data

    tixChainMethod $w ConstructWidget


    frame $w.frame
    set data($w:label) [ label $w.frame.label -text $data(-label) ]
    if {$data(-labelwidth) != ""} {
	$w.frame.label configure -width $data(-labelwidth)
    }
    
    label $w.frame.padlbl -width 1

    frame $w.frame.f -borderwidth 2

    set data(w:x) [ tixControl $w.frame.f.x -label "X" -step 1 -variable $data(-variablex)\
			-selectmode normal -options { entry.width 6 } ]
    set data(w:y) [ tixControl $w.frame.f.y -label "Y" -step 1 -variable $data(-variabley)\
			-selectmode normal -options { entry.width 6 } ]
    set data(w:z) [ tixControl $w.frame.f.z -label "Z" -step 1 -variable $data(-variablez)\
			-selectmode normal -options { entry.width 6 } ]

    radiobutton $w.frame.vox -text "vox" -value vox -variable $data(-variableunits)
    radiobutton $w.frame.mm  -text "mm"  -value mm  -variable $data(-variableunits)
    
    pack $w.frame.f.x $w.frame.f.y $w.frame.f.z -in $w.frame.f -side left -padx 2 -pady 0
 
    pack $w.frame.label $w.frame.padlbl $w.frame.f $w.frame.vox $w.frame.mm -in $w.frame -side left -padx 0 -pady 0

    pack $w.frame -in $w 
}

proc coordinateEdit:update { w } {

    upvar #0 $w data

    [$w subwidget x] update
    [$w subwidget y] update
    [$w subwidget z] update
}

proc mm_to_voxels { X Y Z mask } {

    global FSLDIR

    upvar $X cX
    upvar $Y cY
    upvar $Z cZ


    set vcX [ exec sh -c "echo $cX $cY $cZ | $FSLDIR/bin/std2imgcoord -img $mask -vox - | awk '{print \$1}'" ]    
    set vcY [ exec sh -c "echo $cX $cY $cZ | $FSLDIR/bin/std2imgcoord -img $mask -vox - | awk '{print \$2}'" ] 
    set vcZ [ exec sh -c "echo $cX $cY $cZ | $FSLDIR/bin/std2imgcoord -img $mask -vox - | awk '{print \$3}'" ] 	
    set cX $vcX
    set cY $vcY
    set cZ $vcZ


    #set o1 [ exec sh -c "$FSLDIR/bin/avwval $mask origin1" ]
    #set o2 [ exec sh -c "$FSLDIR/bin/avwval $mask origin2" ]
    #set o3 [ exec sh -c "$FSLDIR/bin/avwval $mask origin3" ]

    #if { $o1 != 0 && $o2 != 0 && $o3 != 0 } {
#	set o1 [ expr $o1 - 0.5 ]
#	set o2 [ expr $o2 - 1.5 ]
#	set o3 [ expr $o3 - 0.5 ]
    #}
#    puts "Using origin $o1 $o2 $o3"

    #set cX [ expr int( ( double($cX) / [ exec sh -c "$FSLDIR/bin/avwval $mask pixdim1" ] ) + $o1 ) ]
    #set cY [ expr int( ( double($cY) / [ exec sh -c "$FSLDIR/bin/avwval $mask pixdim2" ] ) + $o2 ) ]
    #set cZ [ expr int( ( double($cZ) / [ exec sh -c "$FSLDIR/bin/avwval $mask pixdim3" ] ) + $o3 ) ]

#    puts "Setting co-ordinates $cX $cY $cZ"
}

proc fdt:dialog { w tclstartupfile } {

    global tool
    global eddy
    global bedpost
    global registration
    global dtifit
    global probtrack
    global HTMLPATH
    global FSLDIR
    global VERSION
    global INMEDX
    global VARS

    set VARS(history) ""

    set LWIDTH 27
    set PWD [ pwd ]

    if [winfo exists $w] {
        wm deiconify $w
        raise $w
        return
    }

    toplevel $w
    wm title $w "FDT - FMRIB's Diffusion Toolbox $VERSION"
    wm iconname $w "FDT"
    wm iconbitmap $w @$FSLDIR/tcl/fmrib.xbm

    #-------- Stage and Mode Options -------- 

    frame $w.tool
    
    tixOptionMenu $w.tool.menu -variable tool 
    $w.tool.menu add command eddy_current -label "Eddy current correction"
    $w.tool.menu add command bedpost      -label "BEDPOST Estimation of diffusion parameters"
    $w.tool.menu add command registration -label "Registration"
    $w.tool.menu add command probtrack    -label "ProbTrack Probabilistic tracking"
    $w.tool.menu add separator xutilsx
    $w.tool.menu add command dtifit       -label "DTIFit Reconstruct diffusion tensors"

    pack $w.tool.menu -in $w.tool -side left -pady 3 -padx 6 -anchor nw

    #-------- Tool Options... -------- 

    frame $w.opts

    #------- Registration --------

    frame $w.registration

    proc registration_set_directory { w dirname } {
	global registration

	set struct [ file join $dirname struct_brain ]

	if { [ imtest $struct ] } {
	    set registration(struct_image) $struct
	} else {
	    set registration(struct_image) ""
	}
    }

    FSLFileEntry $w.registration.directory \
	-variable registration(directory) \
	-directory $PWD \
	-label "BEDPOST directory:" \
	-labelwidth $LWIDTH \
	-title "Choose directory" \
	-width 35 \
	-command "registration_set_directory $w" \
	-filterhist VARS(history)

    registrationImageSelect $w.registration.struct \
	-variable registration(struct_yn) \
	-filename registration(struct_image) \
	-dof registration(struct_dof) \
	-search registration(struct_search) \
	-costfn registration(struct_costfn) \
	-directory $PWD \
	-label "Main structural image" \
	-labelwidth $LWIDTH \
	-filterhist VARS(history)
    set registration(struct_costfn) mutualinfo
    set registration(struct_dof) 12
    set registration(struct_search) 90

    set registration(standard_yn) 1
    registrationImageSelect $w.registration.standard \
	-variable registration(standard_yn) \
	-filename registration(standard_image) \
	-dof registration(standard_dof) \
	-search registration(standard_search) \
	-costfn registration(standard_costfn) \
	-directory $PWD \
	-label "Standard space" \
	-labelwidth $LWIDTH \
	-filterhist VARS(history)
    set registration(standard_dof) 12
    set registration(standard_search) 90
    set registration(standard_image) [ file join ${FSLDIR} etc standard avg152T1_brain ]

    pack $w.registration.directory $w.registration.struct  $w.registration.standard -in $w.registration -side top -padx 3 -pady 3 -anchor nw
    
    #------- ECC --------

    frame $w.ecc

    proc ecc_update_files { w filename } {
	global eddy

	set eddy(output) [ file join [file dirname $eddy(input)] data ]
    }
    FSLFileEntry $w.ecc.input \
	-variable eddy(input) \
	-directory $PWD \
	-label "Diffusion weighted data:" \
	-labelwidth $LWIDTH \
	-title "Choose diffusion weighted image" \
	-pattern "IMAGE" \
	-width 35 \
	-command "ecc_update_files $w" \
	-filterhist VARS(history)

    FSLFileEntry $w.ecc.output \
	-variable eddy(output) \
	-directory $PWD \
	-label "Corrected output data:" \
	-labelwidth $LWIDTH \
	-title "Choose output image name" \
	-pattern "IMAGE" \
	-width 35 \
	-filterhist VARS(history)
    tixControl $w.ecc.refnum -label "Reference volume" \
	-variable eddy(refnum) -step 1 -min 0 \
	-selectmode immediate -options { entry.width 6 } 
    pack $w.ecc.input $w.ecc.output $w.ecc.refnum -in $w.ecc -side top -padx 3 -pady 3 -expand yes -anchor w

    #------- DTIFit --------

    frame $w.dtifit

    FSLFileEntry $w.dtifit.directory \
	-variable dtifit(directory) \
	-directory $PWD \
	-label "Input directory:" \
	-labelwidth $LWIDTH \
	-title "Choose directory" \
	-command "set_working_directory dtifit(cwd)" \
	-width 35 \
	-filterhist VARS(history)

    proc dtifit_toggle_expert { w } {
	global dtifit

	if { $dtifit(expert_yn) } {
	    pack forget $w.dtifit.directory
	    pack $w.dtifit.expert -in $w.dtifit -after $w.dtifit.expert_yn
	} else {
	    pack forget $w.dtifit.expert
	    pack $w.dtifit.directory -in $w.dtifit -before $w.dtifit.expert_yn
	}
    }

    checkbutton $w.dtifit.expert_yn -text "Specify input files manually" \
	-variable dtifit(expert_yn) -command "dtifit_toggle_expert $w"

    frame $w.dtifit.expert

    proc set_working_directory { cwd filename } {
	upvar $cwd myCWD

	set dirname [file dirname $filename]
	puts "switching from $myCWD to $dirname" 
	set myCWD $dirname
    }

    proc dtifit_update_files { w filename } {
	global dtifit

	set dtifit(output) [ file join [file dirname $dtifit(input)] dti ]
	set_working_directory dtifit(cwd) $dtifit(input)
    }
    
    set dtifit(cwd) $PWD
    FSLFileEntry $w.dtifit.expert.input \
	-variable dtifit(input) \
	-directory $dtifit(cwd) \
	-label "Diffusion weighted data:" \
	-labelwidth $LWIDTH \
	-title "Choose diffusion weighted image" \
	-pattern "IMAGE" \
	-width 35 \
	-command "dtifit_update_files $w" \
	-filterhist VARS(history)
    
    FSLFileEntry $w.dtifit.expert.mask \
	-variable dtifit(mask) \
	-directory $dtifit(cwd) \
	-label "BET binary brain mask:" \
	-labelwidth $LWIDTH \
	-title "Choose BET brain mask file" \
	-pattern "IMAGE" \
	-command "set_working_directory dtifit(cwd)" \
	-width 35 \
	-filterhist VARS(history)
    
    FSLFileEntry $w.dtifit.expert.output \
	-variable dtifit(output) \
	-directory $dtifit(cwd) \
	-label "Output basename:" \
	-labelwidth $LWIDTH \
	-title "Choose output base name" \
	-command "set_working_directory dtifit(cwd)" \
	-width 35 \
	-filterhist VARS(history)
    
    FSLFileEntry $w.dtifit.expert.bvecs \
	-variable dtifit(bvecs) \
	-directory $dtifit(cwd) \
	-label "Gradient directions:" \
	-labelwidth $LWIDTH \
	-title "Choose bvecs file" \
	-command "set_working_directory dtifit(cwd)" \
	-pattern "*" \
	-width 35 \
	-filterhist VARS(history)
    
    FSLFileEntry $w.dtifit.expert.bvals \
	-variable dtifit(bvals) \
	-directory $dtifit(cwd) \
	-label "b values:" \
	-labelwidth $LWIDTH \
	-title "Choose bvals file" \
	-command "set_working_directory dtifit(cwd)" \
	-pattern "*" \
	-width 35 \
	-filterhist VARS(history)
    
    pack $w.dtifit.expert.input $w.dtifit.expert.mask $w.dtifit.expert.output \
	$w.dtifit.expert.bvecs $w.dtifit.expert.bvals \
	-in $w.dtifit.expert -side top -padx 3 -pady 3 -expand yes -anchor w

    pack $w.dtifit.directory $w.dtifit.expert_yn -in $w.dtifit \
	-side top -padx 3 -pady 3 -expand yes -anchor w

    #------- BEDPOST --------

    frame $w.bedpost

    FSLFileEntry $w.bedpost.directory \
	-variable bedpost(directory) \
	-directory $PWD \
	-label "Input directory:" \
	-labelwidth $LWIDTH \
	-title "Choose directory" \
	-width 35 \
	-filterhist VARS(history)

    set bedpost(ecc_yn) 0

    pack $w.bedpost.directory \
	-in $w.bedpost -side top -padx 3 -pady 3 -expand yes -anchor w

    #-------- ...ProbTrack... -------- 

    tixNoteBook $w.probtrack -ipadx 5 -ipady 5

    $w.probtrack add data -label "Data"
    $w.probtrack add options -label "Options"

    #-------- ...Mode specific options... --------

    set dataf [$w.probtrack subwidget data]
    frame $w.data

    tixOptionMenu $w.data.mode -label "Mode: " -variable probtrack(mode)
    [ $w.data.mode subwidget menu ] configure -disabledforeground darkred
    $w.data.mode add command xtitlex -label "Path distribution estimation" -state disabled
    $w.data.mode add command simple  -label "  Single seed voxel"
    $w.data.mode add command all     -label "  Seed mask"
    $w.data.mode add command maska   -label "  Seed mask and waypoint masks"
    $w.data.mode add command masks   -label "  Two masks - symmetric"
    $w.data.mode add separator xseedsx
    $w.data.mode add command seeds  -label "Connectivity-based seed classification"

    if { [ file exists [ file join $FSLDIR bin reord_OM ] ] } {
	$w.data.mode add separator xmatrixx
	$w.data.mode add command mat1   -label "Matrix 1"
	$w.data.mode add command mat2   -label "Matrix 2"
	$w.data.mode add command mskmat -label "Mask matrix"
    }

    set probtrack(xfm) ""
    set probtrack(basename) "merged"
    set probtrack(mask) "nodif_brain_mask"

    proc probtrack_update_files { w filename } {
	global probtrack
	global FSLDIR

	if { ($probtrack(bedpost_dir) != "") && ($probtrack(seed) != "") } {
	    set probtrack(dir) \
		[ file join $probtrack(bedpost_dir) [ file tail [ exec $FSLDIR/bin/remove_ext $probtrack(seed) ] ] ]
	}
    }

    FSLFileEntry $w.data.directory \
	-variable probtrack(bedpost_dir) \
	-directory $PWD \
	-label "BEDPOST directory" \
	-labelwidth $LWIDTH \
	-title "Choose BEDPOST directory" \
	-width 35 \
	-command "probtrack_update_files $w" \
	-filterhist VARS(history)

    tixLabelFrame $w.data.target -label "Target list"
    multiFileSelect $w.data.targets -label "Target file: " -labelwidth $LWIDTH -directory $PWD
    pack $w.data.targets -in [$w.data.target subwidget frame] \
	-side top -anchor w -padx 3 -pady 3

    tixLabelFrame $w.data.seedspace -label "Seed space"
    
    FSLFileEntry $w.data.seed \
	-variable probtrack(seed) \
	-directory $PWD \
	-label "Seed image:" \
	-labelwidth $LWIDTH \
	-title "Choose seed image" \
	-pattern "IMAGE" \
	-width 35 \
	-command "probtrack_update_files $w" \
	-filterhist VARS(history)

    FSLFileEntry $w.data.seed2 \
	-variable probtrack(seed2) \
	-directory $PWD \
	-label "Target image:" \
	-labelwidth $LWIDTH \
	-title "Choose target image" \
	-pattern "IMAGE" \
	-width 35 \
	-command "probtrack_update_files $w" \
	-filterhist VARS(history)

    tixLabelFrame $w.data.waypoint -label "Waypoint masks"
    multiFileSelect $w.data.waypoints -label "Waypoint mask file: " -labelwidth $LWIDTH -directory $PWD
    pack $w.data.waypoints -in [$w.data.waypoint subwidget frame] \
	-side top -anchor w -padx 3 -pady 3

    set probtrack(x) 0
    set probtrack(y) 0
    set probtrack(z) 0
    set probtrack(units) vox

    coordinateEdit $w.data.seedxyz -label "Seed:" -labelwidth $LWIDTH \
	-variablex probtrack(x) \
	-variabley probtrack(y) \
	-variablez probtrack(z) \
	-variableunits probtrack(units)

    proc probtrack_toggle_reference { w } {
	global probtrack

	if { $probtrack(usereference_yn) } { 
	    if { $probtrack(mode) == "simple" } {
		pack $w.data.reference $w.data.xfm -in [$w.data.seedspace subwidget frame] \
		    -side top -after $w.data.usereference_yn -padx 3 -pady 3
	    } else {
		pack $w.data.xfm -in [$w.data.seedspace subwidget frame] \
		    -side top -after $w.data.usereference_yn -padx 3 -pady 3
	    }
	} else { 
	    pack forget $w.data.reference
	    pack forget $w.data.xfm
	}
    }

    checkbutton $w.data.usereference_yn -text "Seed space is not diffusion" \
	-variable probtrack(usereference_yn) \
	-command "probtrack_toggle_reference $w"

    FSLFileEntry $w.data.reference \
	-variable probtrack(reference) \
	-directory $PWD \
	-label "Seed space reference image:" \
	-labelwidth $LWIDTH \
	-title "Choose reference image" \
	-pattern "IMAGE" \
	-width 35 \
	-filterhist VARS(history)

    FSLFileEntry $w.data.xfm \
	-variable probtrack(xfm) \
	-directory $PWD \
	-label "Seed to diff-space transform:" \
	-labelwidth $LWIDTH \
	-title "Select seed-space to DTI-space transformation matrix" \
	-pattern "*.mat" \
	-width 35 \
	-filterhist VARS(history)

    proc probtrack_toggle_exclude { w } {
	global probtrack

	if { $probtrack(exclude_yn) } { 
	    pack $w.data.exclude -in [$w.data.seedspace subwidget frame] \
		-side top -after $w.data.exclude_yn -padx 3 -pady 3
	} else { 
	    pack forget $w.data.exclude
	}
    }

    checkbutton $w.data.exclude_yn -text "Use exclusion mask" -variable probtrack(exclude_yn) \
	-command "probtrack_toggle_exclude $w"

    FSLFileEntry $w.data.exclude \
	-variable probtrack(exclude) \
	-directory $PWD \
	-label "Exclusion mask:" \
	-labelwidth $LWIDTH \
	-title "Select exclusion mask image" \
	-pattern "IMAGE" \
	-width 35 \
	-filterhist VARS(history)

    pack $w.data.seed $w.data.usereference_yn $w.data.exclude_yn -in\
	[$w.data.seedspace subwidget frame] -padx 3 -pady 3 -side top -anchor w

    probtrack_toggle_reference $w
    probtrack_toggle_exclude $w

    FSLFileEntry $w.data.lrmask \
	-variable probtrack(lrmask) \
	-directory $PWD \
	-label "Low resolution mask:" \
	-labelwidth $LWIDTH \
	-title "Choose low resolution mask" \
	-pattern "IMAGE" \
	-width 35 \
	-filterhist VARS(history)

    pack $w.data.lrmask \
	-in [$w.data.seedspace subwidget frame] \
	-side top -padx 3 -pady 3 -anchor w

    tixLabelFrame $w.data.output -label "Ouput"

    FSLFileEntry $w.data.dir \
	-variable probtrack(dir) \
	-directory $PWD \
	-label "Output directory:" \
	-labelwidth $LWIDTH \
	-title "Name the output directory" \
	-width 35 \
	-filterhist VARS(history)

    FSLFileEntry $w.data.out \
	-variable probtrack(output) \
	-directory $PWD \
	-label "Output:" \
	-labelwidth $LWIDTH \
	-title "Choose output file name" \
	-pattern "IMAGE" \
	-width 35 \
	-filterhist VARS(history)

    pack $w.data.out $w.data.dir -in [$w.data.output subwidget frame] \
	-side top -padx 3 -pady 3 -anchor w

    pack $w.data.mode $w.data.directory $w.data.target $w.data.seedspace $w.data.output \
	-in $w.data -side top -padx 3 -pady 3 -anchor nw

    pack $w.data -in $dataf -side top -padx 3 -pady 3 -anchor nw -expand yes -fill both

    #-------- ...Options... --------

    set optionsf [$w.probtrack subwidget options]

    tixLabelFrame $w.options -label "Basic Options"

    checkbutton $w.options.verbose \
	    -text "Verbose" \
	    -variable probtrack(verbose_yn)
    
    set probtrack(nparticles) 5000
    tixControl $w.options.nparticles -label "Number of samples" \
	-variable probtrack(nparticles) -step 100 -min 1 \
	-selectmode immediate -options { entry.width 6 }    
    $w.options.nparticles subwidget label config -width $LWIDTH

    set probtrack(curvature) 0.2
    tixControl $w.options.curvature -label "Curvature threshold" \
	-variable probtrack(curvature) -step 0.01 -min 0.0 -max 1.0 \
	-selectmode immediate -options { entry.width 4 }    
    $w.options.curvature subwidget label config -width $LWIDTH

    set probtrack(loopcheck_yn) 1
    checkbutton $w.options.loopcheck \
	-text "Loopcheck" \
	-variable probtrack(loopcheck_yn)

    collapsible frame $w.advanced -title "Advanced Options"

    set probtrack(nsteps) 2000
    tixControl $w.advanced.nsteps -label "Maximum number of steps" \
	-variable probtrack(nsteps) -step 10 -min 2 \
	-selectmode immediate -options { entry.width 6 }    
    $w.advanced.nsteps subwidget label config -width $LWIDTH

    set probtrack(steplength) 0.5
    tixControl $w.advanced.steplength -label "Step length (mm)" \
	-variable probtrack(steplength) -step 0.1 -min 0 \
	-selectmode immediate -options { entry.width 4 }    
    $w.advanced.steplength subwidget label config -width $LWIDTH

#     set probtrack(usef_yn) 0
#     checkbutton $w.advanced.usef \
# 	-text "Use anisotropy constraints" \
# 	-variable probtrack(usef_yn) \
# 	-command "if { \$probtrack(usef_yn) } { set probtrack(nsteps) 20000 } else { set probtrack(nsteps) 1000 }"

    set probtrack(modeuler_yn) 0
    checkbutton $w.advanced.modeuler \
	-text "Use modified Euler streamlining" \
	-variable probtrack(modeuler_yn)

    pack \
	$w.advanced.modeuler \
	$w.advanced.nsteps \
	$w.advanced.steplength \
	-in $w.advanced.b -side top -pady 3 -padx 6 -anchor nw

    pack \
	$w.options.nparticles \
	$w.options.curvature \
	$w.options.verbose \
	$w.options.loopcheck \
	-in [$w.options subwidget frame] -side top -pady 3 -padx 6 -anchor nw

    pack $w.options $w.advanced -in $optionsf -side top -pady 3 -padx 6 -anchor nw -expand yes -fill both

    #-------- Buttons --------

    frame $w.btns
    frame $w.btns.b -relief raised -borderwidth 1
    
    button $w.apply     -command "fdt:apply $w keep" \
        -text "Go" -width 5
    bind $w.apply <Return> {
        [winfo toplevel %W].apply invoke}
 
    button $w.cancel    -command "fdt:destroy $w" \
        -text "Exit" -width 5
    bind $w.cancel <Return> {
        [winfo toplevel %W].cancel invoke}

    button $w.help -command "FmribWebHelp file: $HTMLPATH/index.html" \
	    -text "Help" -width 5
    bind $w.help <Return> {
	[winfo toplevel %W].help invoke}
 
    pack $w.btns.b -side bottom -fill x -padx 3 -pady 5
    pack $w.apply $w.cancel $w.help -in $w.btns.b \
        -side left -expand yes -padx 3 -pady 10 -fill y
 
    pack $w.tool $w.opts $w.btns -side top -expand yes -fill both
    
    $w.data.mode configure -command "fdt:probtrack_mode $w"
    $w.tool.menu configure -command "fdt:select_tool $w"

    set tool probtrack
    fdt:select_tool $w $tool
    set probtrack(mode) simple
    fdt:probtrack_mode $w $probtrack(mode)

    update idletasks
    if { $tclstartupfile != "" } {
	puts "Reading $tclstartupfile"
	source $tclstartupfile
	fdt:select_tool $w $tool
	fdt:probtrack_mode $w $probtrack(mode)
    }
}

proc fdt:probtrack_mode { w mode } {
    global probtrack

    set seedspacef [$w.data.seedspace subwidget frame]

    pack forget $w.data.target
    pack forget $w.data.lrmask
    pack forget $w.data.usereference_yn
    pack forget $w.data.reference
    pack forget $w.data.seed2
    pack forget $w.data.waypoint
    pack forget $w.data.seedxyz
    pack forget $w.data.out
    pack $w.data.seed $w.data.usereference_yn $w.data.xfm -in $seedspacef -side top \
	-before $w.data.exclude_yn -padx 3 -pady 3 -anchor nw
    [$w.data.seed subwidget label] configure -text "Seed image:"
    pack $w.data.dir -in [$w.data.output subwidget frame] -side top -padx 3 -pady 3 -anchor nw
    switch -- $mode {
  	simple {
 	    pack $w.data.seedxyz -in $seedspacef -side top -padx 3 -pady 3 -before $w.data.usereference_yn
	    pack $w.data.reference -in $seedspacef -side top -padx 3 -pady 3 -before $w.data.xfm -anchor nw
	    pack $w.data.out -in [$w.data.output subwidget frame] -side top -padx 3 -pady 3 -anchor nw
  	    pack forget $w.data.dir
	    pack forget $w.data.seed
  	}
	maska {
	    pack $w.data.waypoint -in $w.data -side top -padx 3 -pady 3 -after $w.data.seedspace -anchor nw
	}
	masks {
	    pack $w.data.seed2 -in $seedspacef -side top -after $w.data.seed
	    [$w.data.seed subwidget label] configure -text "Mask image 1:"
	    [$w.data.seed2 subwidget label] configure -text "Mask image 2:"
	}
  	seeds {
  	    pack $w.data.target -in $w.data -side top -padx 3 -pady 3 -after $w.data.seedspace -anchor nw
  	}
  	mat2 {
  	    pack $w.data.lrmask -in $seedspacef \
 		-side top -padx 3 -pady 3 -after $w.data.xfm -anchor nw
  	}
    }
    probtrack_toggle_reference $w
    probtrack_toggle_exclude $w
}

proc fdt:select_tool { w tool } {
    pack forget $w.ecc
    pack forget $w.probtrack
    pack forget $w.bedpost
    pack forget $w.registration
    pack forget $w.dtifit
    if {$tool == "bedpost"} { 
	pack $w.bedpost -in $w.opts -side top -padx 3 -pady 3 -anchor nw
    } elseif {$tool == "probtrack"}  { 
	pack $w.probtrack -in $w.opts -side top -padx 3 -pady 3 -anchor nw
    } elseif {$tool == "dtifit"}  { 
	pack $w.dtifit -in $w.opts -side top -padx 3 -pady 3 -anchor nw
    } elseif {$tool == "eddy_current"}  { 
	pack $w.ecc -in $w.opts -side top -padx 3 -pady 3 -anchor nw
    } elseif {$tool == "registration"} {
	pack $w.registration -in $w.opts -side top -padx 3 -pady 3 -anchor nw
    }
}
proc fdt_monitor_short { w cmd } {
    global debugging OSFLAVOUR

    puts "$cmd"

    if { $OSFLAVOUR != "cygwin" } {
	set oldcursor [ $w configure -cursor { watch red white } ]

	catch {
	    update idletasks
	    if { ! $debugging } {
		set fd [ open "|$cmd" r ]
#		set fd [ open "|qrsh -V -now n -q short.q $cmd" r ]
		while { ( [ gets $fd line ] >= 0 ) } {
		    update idletasks
		    puts $line
		}
		close $fd
	    }
	} junk

	$w configure -cursor $oldcursor

    } else {
	catch { exec sh -c $cmd } junk
    }

    if { $junk != "" } {
	MxPause "Errors: $junk"
    } 

    puts "Done!"
}

proc fdt_monitor { w cmd } {
    global debugging OSFLAVOUR

    puts "$cmd"

    if { $OSFLAVOUR != "cygwin" } {
	set oldcursor [ $w configure -cursor { watch red white } ]

	catch {
	    update idletasks
	    if { ! $debugging } {
		set fd [ open "|$cmd" r ]
#		set fd [ open "|qrsh -V -now n -q long.q $cmd" r ]
		while { ( [ gets $fd line ] >= 0 ) } {
		    update idletasks
		    puts $line
		}
		close $fd
	    }
	} junk

	$w configure -cursor $oldcursor

    } else {
	catch { exec sh -c $cmd } junk
    }

    if { $junk != "" } {
	MxPause "Errors: $junk"
    } 

    puts "Done!"
}

proc fdt:apply { w dialog } {

    global tool
    global BINPATH
    global FSLDIR

    switch -- $tool {
	eddy_current {
	    global eddy

	    set errorStr ""
	    if { $eddy(input) == "" } { set errorStr "You need to specify the input image! " }
	    if { $eddy(output) == "" } { set errorStr "$errorStr You need to specify an output image!" }
	    if { $errorStr != "" } {
		MxPause $errorStr
		return
	    }

	    #	    check output!=input
	    set canwrite 1
	    if { $eddy(input) == $eddy(output) } {
		set canwrite [ YesNoWidget "Output and input images have the same name. Overwrite input?" Yes No ]
	    }
	    if { $canwrite } {
		fdt_monitor $w "${FSLDIR}/bin/eddy_correct $eddy(input) $eddy(output) $eddy(refnum)"
	    }
	}
	dtifit {
	    global dtifit

	    if { ! $dtifit(expert_yn) } {
		set dtifit(input)  [ file join $dtifit(directory) data ]
		set dtifit(output) [ file join $dtifit(directory) dti ]
		set dtifit(mask)   [ file join $dtifit(directory) nodif_brain_mask ]
		set dtifit(bvecs)  [ file join $dtifit(directory) bvecs ]
		set dtifit(bvals)  [ file join $dtifit(directory) bvals ]
	    }

	    set errorStr ""
	    if { $dtifit(directory) == "" && ! $dtifit(expert_yn) } { set errorStr "You must specify the input directory!" }
	    if { $dtifit(input) == "" } { set errorStr "You need to specify the diffusion weighted data image!" }
	    if { $dtifit(output) == "" } { set errorStr "$errorStr You need to specify the output basename!" }
	    if { $dtifit(mask) == "" } { set errorStr "$errorStr You need to specify a mask image!" }
	    if { $dtifit(bvecs) == "" } { set errorStr "$errorStr Please select a gradient directions file!" }
	    if { $dtifit(bvals) == "" } { set errorStr "$errorStr Please select a b values file!" }
	    if { $errorStr != "" } {
		MxPause $errorStr
		return
	    }

	    set canwrite 1
	    if { [file exists $dtifit(output) ] } {
		set canwrite [ YesNoWidget "Overwrite $dtifit(output)?" Yes No ]
	    }
	    if { $canwrite } {
		fdt_monitor_short $w "${FSLDIR}/bin/dtifit --data=$dtifit(input) --out=$dtifit(output) --mask=$dtifit(mask) --bvecs=$dtifit(bvecs) --bvals=$dtifit(bvals)"
	    }
	}
	bedpost {
	    global bedpost

	    set errorStr ""
	    if { $bedpost(directory) == ""  } { set errorStr "You must specify the bedpost directory!" }
	    if { $errorStr != "" } {
		MxPause $errorStr
		return
	    }

	    set canwrite 1
	    if { [file exists ${bedpost(directory)}.bedpost ] } {
		set canwrite [ YesNoWidget "Overwrite ${bedpost(directory)}.bedpost?" Yes No ]
		if { $canwrite } {
		    puts "rm -rf ${bedpost(directory)}.bedpost"
		    catch { exec rm -rf ${bedpost(directory)}.bedpost } errmsg
		}
	    }
	    if { $canwrite } {
		puts "bedpost $bedpost(directory)"
		fdt_monitor $w "${FSLDIR}/bin/bedpost $bedpost(directory)"
	    }
	}
	probtrack {
	    global probtrack

	    set errorStr ""
	    if { $probtrack(bedpost_dir) == ""  } { set errorStr "You must specify the bedpost directory!" }
	    if { $probtrack(usereference_yn) && $probtrack(xfm) == "" } { set errorStr "$errorStr You must specify the seed-diff transform!" }
	    if { $probtrack(exclude_yn) && $probtrack(exclude) == "" } { set errorStr "$errorStr You must specify the exclusion mask!" }

	    set flags ""
	    if { $probtrack(verbose_yn) == 1 } { set flags "$flags -V 1" }
	    if { $probtrack(loopcheck_yn) == 1 } { set flags "$flags -l" }
#	    if { $probtrack(usef_yn) == 1 } { set flags "$flags -f" }
	    if { $probtrack(modeuler_yn) == 1 } { set flags "$flags --modeuler" }
	    set flags "$flags -c $probtrack(curvature) -S $probtrack(nsteps) --steplength=$probtrack(steplength) -P $probtrack(nparticles)"

	    set tn [open "| $BINPATH/tmpnam"]
	    gets $tn filebase
	    close $tn
	    set logfile "${filebase}_log.tcl"
	    set log [open "$logfile" w]
	    puts $log "set tool $tool"
	    set copylog ""
	    set ssopts ""
	    if { $probtrack(usereference_yn) } {
		set ssopts "--xfm=$probtrack(xfm)"
	    }
	    if { $probtrack(exclude_yn) == 1 } {
		set ssopts "$ssopts --rubbish=$probtrack(exclude)"
		puts $log "set probtrack(exclude) $probtrack(exclude)"
	    }
	    set basics "--forcedir -s $probtrack(bedpost_dir)/merged -m $probtrack(bedpost_dir)/nodif_brain_mask"	    

    	    foreach entry {bedpost_dir xfm mode exclude_yn usereference_yn verbose_yn loopcheck_yn modeuler_yn \
			       curvature nsteps steplength nparticles} {
		puts $log "set probtrack($entry) $probtrack($entry)"
	    }

	    switch -- $probtrack(mode) {
		simple {
		    if { $probtrack(output) == ""  } { set errorStr "$errorStr You must specify the output basename!" }
		    if { $probtrack(usereference_yn) && $probtrack(reference) == "" } { 
			set errorStr "$errorStr You must specify a seed space reference image!" 
		    }
		    if { $probtrack(usereference_yn) == 1 } {
			set ssopts "--xfm=$probtrack(xfm) --seedref=$probtrack(reference)"
			puts $log "set probtrack(reference) $probtrack(reference)"
			puts $log "set probtrack(xfm) $probtrack(xfm)"
		    } else {
			set ssopts ""
		    }
		    if { $probtrack(exclude_yn) == 1 } {
			set ssopts "$ssopts --rubbish=$probtrack(exclude)"
		    }
		    set fd [ open "${filebase}_coordinates.txt" w ]
		    $w.data.seedxyz update
		    if { $probtrack(units) == "mm" } {
			set x $probtrack(x)
			set y $probtrack(y)
			set z $probtrack(z)
			if { $probtrack(usereference_yn) && $probtrack(reference) != "" } {
			    mm_to_voxels x y z $probtrack(reference)
			} else {
			    mm_to_voxels x y z [ file join $probtrack(bedpost_dir) nodif_brain_mask ]
			}			    
			puts $fd "$x $y $z"
			puts "$probtrack(x) $probtrack(y) $probtrack(z) (mm) -> $x $y $z (voxels)"
		    } else {
			puts $fd "$probtrack(x) $probtrack(y) $probtrack(z)"
		    }
		    close $fd
 		    puts $log "set probtrack(x) $probtrack(x)"
		    puts $log "set probtrack(y) $probtrack(y)"
		    puts $log "set probtrack(z) $probtrack(z)"
		    puts $log "set probtrack(units) $probtrack(units)"
		    puts $log "set probtrack(output) $probtrack(output)"

		    if { $errorStr != "" } {
			MxPause $errorStr
			return
		    }

		    set canwrite 1
		    if { [ file exists $probtrack(output) ] } {
			set canwrite [  YesNoWidget "Overwrite $probtrack(output)?" Yes No ]
			if { $canwrite } {
			    puts "rm -rf $probtrack(output)"
			    exec rm -rf $probtrack(output)
			}
		    }
		    if { $canwrite } {
			set copylog "fdt.log"
			fdt_monitor_short $w "$FSLDIR/bin/probtrack --mode=simple -x ${filebase}_coordinates.txt $basics $ssopts $flags -o $probtrack(output)"
		    }
		    puts "rm ${filebase}_coordinates.txt"
		    exec rm ${filebase}_coordinates.txt
		}
		all {
		    if { $probtrack(seed) == ""  } { set errorStr "$errorStr You must specify a seed image!" }
		    if { $probtrack(dir) == ""  } { set errorStr "$errorStr You must specify the output directory!" }
		    if { $errorStr != "" } {
			MxPause $errorStr
			return
		    }
		    puts $log "set probtrack(seed) $probtrack(seed)"
		    puts $log "set probtrack(dir) $probtrack(dir)"

		    set canwrite 1
		    if { [ file exists $probtrack(dir) ] } {
			set canwrite [  YesNoWidget "Overwrite $probtrack(dir)?" Yes No ]
			if { $canwrite } {
			    puts "rm -rf $probtrack(dir)"
			    exec rm -rf $probtrack(dir)
			}
		    }
		    if { $canwrite } {
			puts "mkdir -p $probtrack(dir)"
			exec mkdir -p $probtrack(dir)
			set copylog "$probtrack(dir)/fdt.log"
			fdt_monitor $w "$FSLDIR/bin/probtrack --mode=seedmask -x $probtrack(seed) $basics $ssopts $flags -o fdt_paths --dir=$probtrack(dir)"
		    }
		}
		masks {
		    if { $probtrack(seed2) == "" } { set errorStr "$errorStr You must select a target image!" }
		    if { $probtrack(dir) == ""  } { set errorStr "$errorStr You must specify the output directory!" }
		    if { $errorStr != "" } {
			MxPause $errorStr
			return
		    }
		    puts $log "set probtrack(seed) $probtrack(seed)"
		    puts $log "set probtrack(seed2) $probtrack(seed2)"
		    puts $log "set probtrack(dir) $probtrack(dir)"

		    set canwrite 1
		    if { [ file exists $probtrack(dir) ] } {
			set canwrite [  YesNoWidget "Overwrite $probtrack(dir)?" Yes No ]
			if { $canwrite } {
			    puts "rm -rf $probtrack(dir)"
			    exec rm -rf $probtrack(dir)
			}
		    }
		    if { $canwrite } {
			puts "mkdir -p $probtrack(dir)"
			exec mkdir -p $probtrack(dir)
			set copylog "$probtrack(dir)/fdt.log"
			fdt_monitor $w "$FSLDIR/bin/probtrack --mode=twomasks_symm --seed=$probtrack(seed) --mask2=$probtrack(seed2) $basics $ssopts $flags -o fdt_paths --dir=$probtrack(dir)"
		    }
		}
		maska {
		    if { $probtrack(seed) == ""  } { set errorStr "$errorStr You must specify a seed image!" }
		    if { $probtrack(dir) == ""  } { set errorStr "$errorStr You must specify the output directory!" }
		    if { $errorStr != "" } {
			MxPause $errorStr
			return
		    }
		    puts $log "set probtrack(seed) $probtrack(seed)"
		    puts $log "set probtrack(dir) $probtrack(dir)"

		    set canwrite 1
		    if { [ file exists $probtrack(dir) ] } {
			set canwrite [  YesNoWidget "Overwrite $probtrack(dir)?" Yes No ]
			if { $canwrite } {
			    puts "rm -rf $probtrack(dir)"
			    exec rm -rf $probtrack(dir)
			}
		    }
		    if { $canwrite } {
			puts "mkdir -p $probtrack(dir)"
			exec mkdir -p $probtrack(dir)
			puts $log "$w.data.watpoints load \"$probtrack(dir)/waypoints.txt\""
			$w.data.waypoints save "$probtrack(dir)/waypoints.txt"
			set copylog "$probtrack(dir)/fdt.log"
			fdt_monitor $w "$FSLDIR/bin/probtrack --mode=waypoints --seed=$probtrack(seed) --mask2=${probtrack(dir)}/waypoints.txt $basics $ssopts $flags -o fdt_paths --dir=$probtrack(dir)"
		    }
		}
		seeds {
		    if { $probtrack(seed) == ""  } { set errorStr "$errorStr You must specify a seed image!" }
		    if { $probtrack(dir) == ""  } { set errorStr "$errorStr You must specify the output directory!" }
		    if { $errorStr != "" } {
			MxPause $errorStr
			return
		    }
		    puts $log "set probtrack(seed) $probtrack(seed)"
		    puts $log "set probtrack(dir) $probtrack(dir)"


		    set canwrite 1
		    if { [ file exists $probtrack(dir) ] } {
			set canwrite [  YesNoWidget "Overwrite $probtrack(dir)?" Yes No ]
			if { $canwrite } {
			    puts "rm -rf $probtrack(dir)"
			    exec rm -rf $probtrack(dir)
			}
		    }
		    if { $canwrite } {
			puts "mkdir -p $probtrack(dir)"
			exec mkdir -p $probtrack(dir)
			puts $log "$w.data.targets load \"$probtrack(dir)/targets.txt\""
			$w.data.targets save "$probtrack(dir)/targets.txt"
			set copylog "$probtrack(dir)/fdt.log"
			fdt_monitor $w "$FSLDIR/bin/probtrack --mode=seeds_to_targets -x $probtrack(seed) $basics $ssopts $flags --targetmasks=${probtrack(dir)}/targets.txt --dir=$probtrack(dir)"
		    }
		}
		mat1 {
		    if { $probtrack(seed) == ""  } { set errorStr "$errorStr You must specify a seed image!" }
		    if { $probtrack(dir) == ""  } { set errorStr "$errorStr You must specify the output directory!" }
		    if { $errorStr != "" } {
			MxPause $errorStr
			return
		    }
		    puts $log "set probtrack(seed) $probtrack(seed)"
		    puts $log "set probtrack(dir) $probtrack(dir)"

		    set canwrite 1
		    if { [ file exists $probtrack(dir) ] } {
			set canwrite [  YesNoWidget "Overwrite $probtrack(dir)?" Yes No ]
			if { $canwrite } {
			    puts "rm -rf $probtrack(dir)"
			    exec rm -rf $probtrack(dir)
			}
		    }
		    if { $canwrite } {
			puts "mkdir -p $probtrack(dir)"
			exec mkdir -p $probtrack(dir)
			set copylog "$probtrack(dir)/fdt.log"
			fdt_monitor $w "$FSLDIR/bin/probtrack --mode=matrix1 -x $probtrack(seed) $basics $ssopts $flags -o fdt_matrix --dir=$probtrack(dir)"
		    }
		}
		mat2 {
		    if { $probtrack(seed) == ""  } { set errorStr "$errorStr You must specify a seed image!" }
		    if { $probtrack(dir) == ""  } { set errorStr "$errorStr You must specify the output directory!" }
		    if { $probtrack(lrmask) == "" } { set errorStr "$errorStr You must specify the low resolution mask!" }
		    if { $errorStr != "" } {
			MxPause $errorStr
			return
		    }
		    puts $log "set probtrack(seed) $probtrack(seed)"
		    puts $log "set probtrack(dir) $probtrack(dir)"
		    puts $log "set probtrack(lrmask) $probtrack(lrmask)"

		    set canwrite 1
		    if { [ file exists $probtrack(dir) ] } {
			set canwrite [  YesNoWidget "Overwrite $probtrack(dir)?" Yes No ]
			if { $canwrite } {
			    puts "rm -rf $probtrack(dir)"
			    exec rm -rf $probtrack(dir)
			}
		    }
		    if { $canwrite } {
			puts "mkdir -p $probtrack(dir)"
			exec mkdir -p $probtrack(dir)
			set copylog "$probtrack(dir)/fdt.log"
			fdt_monitor $w "$FSLDIR/bin/probtrack --mode=matrix2 -x $probtrack(seed) $basics $ssopts --lrmask=$probtrack(lrmask) $flags -o fdt_matrix --dir=$probtrack(dir)"
		    }
		}
		mskmat {
		    if { $probtrack(seed) == ""  } { set errorStr "$errorStr You must specify a seed image!" }
		    if { $probtrack(dir) == ""  } { set errorStr "$errorStr You must specify the output directory!" }
		    if { $errorStr != "" } {
			MxPause $errorStr
			return
		    }
		    puts $log "set probtrack(seed) $probtrack(seed)"
		    puts $log "set probtrack(dir) $probtrack(dir)"

		    set canwrite 1
		    if { [ file exists $probtrack(dir) ] } {
			set canwrite [  YesNoWidget "Overwrite $probtrack(dir)?" Yes No ]
			if { $canwrite } {
			    puts "rm -rf $probtrack(dir)"
			    exec rm -rf $probtrack(dir)
			}
		    }
		    if { $canwrite } {
			puts "mkdir -p $probtrack(dir)"
			exec mkdir -p $probtrack(dir)
			set copylog "$probtrack(dir)/fdt.log"
			fdt_monitor $w "$FSLDIR/bin/probtrack --mode=maskmatrix -x $probtrack(seed) $basics $ssopts $flags -o fdt_matrix --dir=$probtrack(dir)"
		    }
		}
	    }
	    close $log
	    if { $copylog != "" } {
		puts "mv $logfile $copylog"
		exec mv $logfile $copylog
	    } else {
		puts "rm $logfile"
		exec rm $logfile
	    }
	}
	registration {
	    global registration

	    set errorStr ""
	    if { $registration(directory) == ""  } { set errorStr "You must specify the bedpost directory!" }
	    if { $registration(struct_yn) && $registration(struct_image) == ""  } { set errorStr "$errorStr You must specify the structural image!" }
	    if { $registration(standard_yn) && $registration(standard_image) == ""  } { set errorStr "$errorStr You must specify the standard image!" }
	    if { $errorStr != "" } {
		MxPause $errorStr
		return
	    }

	    exec mkdir -p [ file join $registration(directory) xfms ]
	    set eyefd [ open [ file join $registration(directory) xfms eye.mat ] w ]
	    puts $eyefd "1 0 0 0"
	    puts $eyefd "0 1 0 0"
	    puts $eyefd "0 0 1 0"
	    puts $eyefd "0 0 0 1"
	    close $eyefd

	    set diff2str   [ file join $registration(directory) xfms diff2str.mat ]
	    set str2diff   [ file join $registration(directory) xfms str2diff.mat ]
	    set str2stand  [ file join $registration(directory) xfms str2standard.mat ]
	    set stand2str  [ file join $registration(directory) xfms standard2str.mat ]
	    set diff2stand [ file join $registration(directory) xfms diff2standard.mat ]
	    set stand2diff [ file join $registration(directory) xfms standard2diff.mat ]
	    set diff       [ file join $registration(directory) nodif_brain ]
	    if { $registration(struct_yn) } {
		set searchrx  "-searchrx -$registration(struct_search) $registration(struct_search)"
		set searchry  "-searchry -$registration(struct_search) $registration(struct_search)"
		set searchrz  "-searchrz -$registration(struct_search) $registration(struct_search)"
		set options   "$searchrx $searchry $searchrz -dof $registration(struct_dof)"
		fdt_monitor $w "${FSLDIR}/bin/flirt -in $diff -ref $registration(struct_image) -omat $diff2str $options -cost $registration(struct_costfn)"
		fdt_monitor $w "${FSLDIR}/bin/convert_xfm -omat $str2diff -inverse $diff2str"
		if { $registration(standard_yn) } {
		    set searchrx  "-searchrx -$registration(standard_search) $registration(standard_search)"
		    set searchry  "-searchry -$registration(standard_search) $registration(standard_search)"
		    set searchrz  "-searchrz -$registration(standard_search) $registration(standard_search)"
		    set options   "$searchrx $searchry $searchrz -dof $registration(standard_dof)"
		    fdt_monitor $w "${FSLDIR}/bin/flirt -in $registration(struct_image) -ref $registration(standard_image) -omat $str2stand $options -cost $registration(standard_costfn)"
		    fdt_monitor $w "${FSLDIR}/bin/convert_xfm -omat $stand2str -inverse $str2stand"
		    fdt_monitor $w "${FSLDIR}/bin/convert_xfm -omat $diff2stand -concat $str2stand $diff2str"
		    fdt_monitor $w "${FSLDIR}/bin/convert_xfm -omat $stand2diff -inverse $diff2stand"
		}
	    } elseif { $registration(standard_yn) } {
		set searchrx  "-searchrx -$registration(standard_search) $registration(standard_search)"
		set searchry  "-searchry -$registration(standard_search) $registration(standard_search)"
		set searchrz  "-searchrz -$registration(standard_search) $registration(standard_search)"
		set options   "$searchrx $searchry $searchrz -dof $registration(standard_dof)"
		fdt_monitor $w "${FSLDIR}/bin/flirt -in $diff -ref $registration(standard_image) -omat $diff2stand $options"
		fdt_monitor $w "${FSLDIR}/bin/convert_xfm -omat $stand2diff -inverse $diff2stand"
	    }
	    puts "Done!"
	    # Fudge to make the logic work
	    set canwrite 1
	}
    }

    if { $canwrite } { 
	MxPause "  Done!  "
	update idletasks
    }

    if {$dialog == "destroy"} {
        fdt:destroy $w
    }
}


proc fdt:destroy { w } {
    destroy $w
}    

set debugging 0

while {[llength $argv] > 0 } {
    set flag [lindex $argv 0]
    switch -- $flag {
	"-debugging" {
	    set debugging 1
	    set argv [lrange $argv 1 end]
	    puts "Debug mode!"
	}
	default { break }
    }
}


wm withdraw .
if { [ info exists env(MRDATADIR) ] } {
    set MRDATADIR $env(MRDATADIR)
} else {
    set MRDATADIR ~/MRdata
}

fdt:dialog .fdt $argv
tkwait window .fdt
