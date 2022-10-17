#

# ApplyXFM - the GUI for applying an xfm
#
# Mark Jenkinson and Stephen Smith, FMRIB Image Analysis Group
#
# Copyright (C) 2001-2006 University of Oxford
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



proc applyxfm { w } {


    global entries FSLDIR PWD

    # ---- Set up Frames ----
    toplevel $w
    wm title $w "ApplyXFM"
    wm iconname $w "ApplyXFM"
    wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm
    tixBalloon    $w.bhelp
    frame $w.f


    tixLabelFrame $w.f.input -label "Input"
    set lfinput [ $w.f.input subwidget frame ]

# Input image

set entries($w,invol) ""

FSLFileEntry $w.f.invol \
	-variable entries($w,invol) \
	-pattern "IMAGE" \
	-directory $PWD \
	-label "Input Volume (3D or 4D)  " \
	-labelwidth 29 \
	-title "Select" \
	-width 40 \
	-filterhist VARS(history)

# Transform matrix

set entries($w,transmat) ""

FSLFileEntry $w.f.xfm \
	-variable entries($w,transmat) \
	-pattern "*.mat" \
	-directory $PWD \
	-label "Transformation Matrix   " \
		-labelwidth 29 \
		-title "Select" \
		-width 40 \
		-filterhist VARS(history)

# Identity and Inverse button
frame $w.f.xfminv

label $w.f.xfminv.idlabel -text "Use Identity transformation       "
set entries($w,idxfm) 0
checkbutton $w.f.xfminv.idbutton -variable entries($w,idxfm) -command "applyxfm:updateinput $w $lfinput "

label $w.f.xfminv.label -text "Apply inverse transformation"
set entries($w,invxfm) 0
checkbutton $w.f.xfminv.button -variable entries($w,invxfm)

  # pack
    pack $w.f.xfminv.idbutton $w.f.xfminv.idlabel $w.f.xfminv.button $w.f.xfminv.label -in $w.f.xfminv -side left -padx 3 -pady 3
    pack $w.f.xfminv $w.f.xfm $w.f.invol -in $lfinput -side top -anchor w -pady 3 -padx 5

# Reference volume

set entries($w,refvol) ""

    tixLabelFrame $w.f.outsize -label "Output Size"
    set lfoutsize [ $w.f.outsize subwidget frame ]

FSLFileEntry $w.f.outsize.refvol \
	-variable entries($w,refvol) \
	-pattern "IMAGE" \
	-directory $PWD \
	-label "Reference Volume   " \
	-labelwidth 29 \
	-title "Select" \
	-width 40 \
	-filterhist VARS(history)


 #output volume size

set entries($w,refsize) 0

set entries($w,outsize_nx) 64
set entries($w,outsize_ny) 64
set entries($w,outsize_nz) 25
set entries($w,outsize_dx) 4
set entries($w,outsize_dy) 4
set entries($w,outsize_dz) 6

    tixOptionMenu $w.f.outsize.volname -label "Base on: " -variable entries($w,refsize)  -command "applyxfm:updateoutsize $w $lfoutsize"
    $w.f.outsize.volname add command 0 -label "Existing Volume"
    $w.f.outsize.volname add command 1 -label "Voxel Dimensions"

    frame $w.f.outsize.n
    frame $w.f.outsize.d

    label $w.f.outsize.n.lab -text "Number of Voxels: " -width 18

    tixControl $w.f.outsize.n.x -label " X " \
	    -variable entries($w,outsize_nx) -step 1 -min 1 -max 10000 \
	    -selectmode immediate
    
    tixControl $w.f.outsize.n.y -label " Y " \
	    -variable entries($w,outsize_ny) -step 1 -min 1 -max 10000 \
	    -selectmode immediate
    
    tixControl $w.f.outsize.n.z -label " Z " \
	    -variable entries($w,outsize_nz) -step 1 -min 1 -max 10000 \
	    -selectmode immediate
    
    label $w.f.outsize.d.lab -text "Voxel Size (mm): " -width 18

    tixControl $w.f.outsize.d.x -label " X " \
	    -variable entries($w,outsize_dx) -step 0.1 -min 0 -max 10000 \
	    -selectmode immediate
    
    tixControl $w.f.outsize.d.y -label " Y " \
	    -variable entries($w,outsize_dy) -step 0.1 -min 0 -max 10000 \
	    -selectmode immediate
    
    tixControl $w.f.outsize.d.z -label " Z " \
	    -variable entries($w,outsize_dz) -step 0.1 -min 0 -max 10000 \
	    -selectmode immediate
    
    pack $w.f.outsize.n.lab $w.f.outsize.n.x $w.f.outsize.n.y $w.f.outsize.n.z -in $w.f.outsize.n -side left -anchor w -padx 3 -pady 3
    pack $w.f.outsize.d.lab $w.f.outsize.d.x $w.f.outsize.d.y $w.f.outsize.d.z -in $w.f.outsize.d -side left -anchor w -padx 3 -pady 3
    pack $w.f.outsize.volname -in $lfoutsize -side top -anchor w -padx 3 -pady 3
#    pack $w.f.outsize.n $w.f.outsize.d -in  $lfoutsize -side top -anchor w -padx 3 -pady 3
    pack  $w.f.outsize.refvol -in  $lfoutsize -side top -anchor w -padx 3 -pady 3



# output volume

set entries($w,outvol) ""

    tixLabelFrame $w.f.output -label "Output"
    set lfoutput [ $w.f.output subwidget frame ]

FSLFileEntry $w.f.outvol \
	-variable entries($w,outvol) \
	-pattern "IMAGE" \
	-directory $PWD \
	-label "Output Volume   " \
	-labelwidth 29 \
		-title "Select" \
		-width 40 \
		-filterhist VARS(history)


    pack $w.f.outvol -in $lfoutput -side top -anchor w -pady 3 -padx 5

    pack $w.f.input $w.f.outsize $w.f.output -in $w.f -side top -anchor w -pady 0 -padx 5


    # advanced options

    # ---- Optional stuff ----

    collapsible frame $w.f.opts -title "Advanced Options"    

    tixNoteBook $w.nb -ipadx 5 -ipady 5

    $w.nb add interp -label "Interpolation Method"
    $w.nb add misc -label "Miscellaneous"
    
    # Interpolation

    set interplf [$w.nb subwidget interp]

    radiobutton $w.trilinear -text "Tri-Linear" \
	    -variable entries($w,interp) -value trilinear -anchor w -command "applyxfm:updateinterp $w $interplf"
    radiobutton $w.nearestneighbour -text "Nearest Neighbour" \
	    -variable entries($w,interp) -value nearestneighbour -anchor w -command "applyxfm:updateinterp $w $interplf"
    radiobutton $w.sinc -text "Sinc" \
	    -variable entries($w,interp) -value sinc -anchor w -command "applyxfm:updateinterp $w $interplf"

    tixControl $w.sincwidth -label " Width of Sinc Window (full width - voxels)" \
	    -variable entries($w,sincwidth) -step 1 -min 1 -max 5000 -selectmode immediate
    set entries($w,sincwidth) 7

    frame $w.swinopt
    label $w.swinbanner -text "Sinc Window Options"
    radiobutton $w.rectangular -text "Rectangular" \
	    -variable entries($w,sincwindow) -value rectangular -anchor w
    radiobutton $w.hanning -text "Hanning" \
	    -variable entries($w,sincwindow) -value hanning -anchor w
    radiobutton $w.blackman -text "Blackman" \
	    -variable entries($w,sincwindow) -value blackman -anchor w
    set entries($w,sincwindow) hanning
    
    # ---- pack ----
    pack $w.trilinear -in $interplf -side top -anchor w -padx 3
    pack $w.nearestneighbour $w.sinc -in $interplf -side top -anchor w -padx 3
    set entries($w,interp) trilinear

    pack $w.swinbanner -in $w.swinopt -side top -anchor w -padx 3
    pack $w.rectangular $w.hanning $w.blackman -in $w.swinopt -side left -anchor w -padx 3


    # Misc

    set misclf [$w.nb subwidget misc]

    set entries($w,datatype) 0
    set entries($w,paddingsize) 0.0

    tixControl $w.paddingsize -label " Extra Padding for Output FOV (in voxels) " \
	    -variable entries($w,paddingsize) -step 0.5 -min 0 -max 5000 -selectmode immediate

    checkbutton $w.datatype  -text " Force Output Datatype to be Floating Point " -variable entries($w,datatype)

    # ---- pack ----

    pack $w.paddingsize $w.datatype -in $misclf -side top -anchor w -padx 3

    set entries($w,datatype) 0
    set entries($w,paddingsize) 0.0


    # ---- pack ----

    frame $w.f.advopts
    pack $w.nb -in $w.f.advopts -side top
    pack $w.f.advopts -in $w.f.opts.b -side left -padx 8 -pady 6 -expand yes -fill both
    pack $w.f.opts -in $w.f -side left -padx 5 -pady 5



    # ---- Button Frame ----

    frame $w.btns
    frame $w.btns.b -relief raised -borderwidth 1
    
    button $w.go     -command "ApplyXFM:go $w" \
	    -text "OK" -width 5

    button $w.apply     -command "ApplyXFM:apply $w" \
	    -text "Apply" -width 5

    button $w.cancel    -command "destroy $w" \
	    -text "Exit" -width 5

    button $w.help -command "FmribWebHelp file: ${FSLDIR}/doc/flirt/overview.html" \
	    -text "Help" -width 5

    pack $w.btns.b -side bottom -fill x
    pack $w.go $w.apply $w.cancel $w.help -in $w.btns.b \
	    -side left -expand yes -padx 3 -pady 10 -fill y
    
    pack $w.f $w.btns -expand yes -fill both


}



proc ApplyXFM:go { w } {
    global entries

    set status [ ApplyXFM:apply $w ]
    destroy $w
}


proc ApplyXFM:apply { w } {
    global entries

    set status [ applyxfm:proc $entries($w,invol) $entries($w,transmat) $entries($w,outvol) $entries($w,refvol) $entries($w,invxfm) $entries($w,refsize) $entries($w,outsize_nx) $entries($w,outsize_ny) $entries($w,outsize_nz) $entries($w,outsize_dx) $entries($w,outsize_dy) $entries($w,outsize_dz) $entries($w,interp) $entries($w,sincwidth) $entries($w,sincwindow) $entries($w,datatype) $entries($w,paddingsize) $entries($w,idxfm) ]

    update idletasks
    puts "Done"
}

proc applyxfm:updateinput { w lfinput } {
    global entries
    
    if { $entries($w,idxfm) == 1 } {
	pack forget $w.f.xfm
	pack forget $w.f.xfminv.button
	pack forget $w.f.xfminv.label
	pack $w.f.xfminv.idbutton $w.f.xfminv.idlabel -in $w.f.xfminv -side left -padx 3 -pady 3
	pack $w.f.xfminv $w.f.invol -in $lfinput -side top -anchor w -pady 3 -padx 5
    } else {
	pack $w.f.xfminv.idbutton $w.f.xfminv.idlabel $w.f.xfminv.button $w.f.xfminv.label -in $w.f.xfminv -side left -padx 3 -pady 3
	pack $w.f.xfminv $w.f.xfm $w.f.invol -in $lfinput -side top -anchor w -pady 3 -padx 5
    }

}


proc applyxfm:updateoutsize { w lfoutsize dummy } {
    global entries
    
    if { $entries($w,refsize) == 1 } {
	pack forget $w.f.outsize.refvol
	pack $w.f.outsize.n $w.f.outsize.d -in  $lfoutsize -side top -anchor w -padx 3 -pady 3
    } else {
	pack forget $w.f.outsize.n $w.f.outsize.d
	pack $w.f.outsize.refvol -in  $lfoutsize -side top -anchor w -padx 3 -pady 3
    }

}


proc applyxfm:updateinterp { w interplf } {
    global entries

    if { [ string match $entries($w,interp) "sinc" ] == 1 } {
	pack $w.swinopt -in $interplf -side top -anchor w -padx 40
	pack $w.sincwidth -in $interplf -side top -anchor w -padx 40
    } else {
	pack forget $w.swinopt
	pack forget $w.sincwidth
    }
}


proc applyxfm:proc { invol transmat outvol refvol invxfm refsize nx ny nz dx dy dz interp sincwidth sincwindow datatype paddingsize idxfm } {

    global FSLDIR
    
    set tmpfiles ""

    if { $idxfm == 1 } {
	set thetransmat $FSLDIR/etc/flirtsch/ident.mat
    } else {
	set thetransmat $transmat
    }

    if { $idxfm == 0 } {
	if { $invxfm == 1 } {
	    set thetransmat ${transmat}_tmp
	    set tmpfiles "$tmpfiles $thetransmat"
	    set invcommand "${FSLDIR}/bin/convert_xfm -omat $thetransmat -inverse $transmat"
	    puts $invcommand
	    catch { exec sh -c $invcommand } errmsg
	    puts $errmsg
	}
    }

    set flirtcommand "${FSLDIR}/bin/flirt -in $invol -applyxfm -init $thetransmat -out $outvol -paddingsize $paddingsize -interp $interp"    

    if { $interp == "sinc" } {
	set flirtcommand "$flirtcommand -sincwidth $sincwidth -sincwindow $sincwindow"
    }

    if { $datatype == 1 } { 
	set flirtcommand "$flirtcommand -datatype float"
    }

    if { $refsize == 1 } {
	set dtype [ exec sh -c "${FSLDIR}/bin/avwval $invol datatype" ]
	set tmpnm ${outvol}_tmp
	set tmpfiles "$tmpfiles $tmpnm"
	set flirtcommand "${FSLDIR}/bin/avwcreatehd $nx $ny $nz 1 $dx $dy $dz 1 0 0 0 $dtype ${tmpnm}.nii.gz ; $flirtcommand -ref ${tmpnm}"
    } else {
	set flirtcommand "$flirtcommand -ref $refvol"
    }

    puts $flirtcommand
    catch { exec sh -c $flirtcommand } errmsg
    puts $errmsg

    # clean up
    # puts "rm $tmpfiles"
    catch { exec sh -c "rm -f $tmpfiles" }

    if { [ imtest $outvol ] == 0 } {
	puts "No output saved!"
	return 4
    }

    return 0
}




wm withdraw .
applyxfm .rename
tkwait window .rename

