#

# InvertXFM - the GUI for convert_xfm
#
# Mark Jenkinson and Stephen Smith, FMRIB Image Analysis Group
#
# Copyright (C) 2001 University of Oxford
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

#{{{ setup

source [ file dirname [ info script ] ]/fslstart.tcl

set VARS(history) {}

#}}}
#{{{ invertxfm

proc invertxfm { w } {

    #{{{ setup main window etc

    global entries FSLDIR PWD

    # ---- Set up Frames ----
    toplevel $w
    wm title $w "InvertXFM"
    wm iconname $w "InvertXFM"
    wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm
    tixBalloon    $w.bhelp
    frame $w.f

#}}}

    tixLabelFrame $w.f.input -label "Input"
    set lfinput [ $w.f.input subwidget frame ]
    #{{{ input transform

set entries($w,1) ""

FSLFileEntry $w.f.xfm \
	-variable entries($w,1) \
	-pattern "*.mat" \
	-directory $PWD \
	-label "Transformation File for A to B   " \
		-labelwidth 29 \
		-title "Select" \
		-width 40 \
		-filterhist VARS(history)

#}}}
    pack $w.f.xfm -in $lfinput -side top -anchor w -pady 3 -padx 5

    tixLabelFrame $w.f.output -label "Output"
    set lfoutput [ $w.f.output subwidget frame ]
    #{{{ output filename

set entries($w,4) ""

FSLFileEntry $w.f.oxfm \
	-variable entries($w,4) \
	-pattern "*.mat" \
	-directory $PWD \
	-label "Save Inverse Transform (B to A)" \
	-labelwidth 29 \
		-title "Select" \
		-width 40 \
		-filterhist VARS(history)

#}}}
    pack $w.f.oxfm -in $lfoutput -side top -anchor w -pady 3 -padx 5

    pack $w.f.input $w.f.output -in $w.f -side top -anchor w -pady 0 -padx 5

    #{{{ buttons

    # ---- Button Frame ----

    frame $w.btns
    frame $w.btns.b -relief raised -borderwidth 1
    
    button $w.go     -command "InvertXFM:go $w" \
	    -text "OK" -width 5

    button $w.apply     -command "InvertXFM:apply $w" \
	    -text "Apply" -width 5

    button $w.cancel    -command "destroy $w" \
	    -text "Exit" -width 5

    button $w.help -command "FmribWebHelp file: ${FSLDIR}/doc/flirt/overview.html" \
	    -text "Help" -width 5

    pack $w.btns.b -side bottom -fill x
    pack $w.go $w.apply $w.cancel $w.help -in $w.btns.b \
	    -side left -expand yes -padx 3 -pady 10 -fill y
    
    pack $w.f $w.btns -expand yes -fill both

#}}}
}

#}}}
#{{{ InvertXFM:apply

proc InvertXFM:go { w } {
    global entries

    catch { InvertXFM:apply $w }
    destroy $w
}

proc InvertXFM:apply { w } {
    global entries

    catch { invertxfm:proc $entries($w,1) $entries($w,4) }

    update idletasks
    puts "Done"
}

#}}}
#{{{ invertxfm:proc

proc invertxfm:proc { transAB invxfmfilename } {

    global FSLDIR

    set thecommand "${FSLDIR}/bin/convert_xfm -omat $invxfmfilename -inverse $transAB"
    puts $thecommand
    catch { exec sh -c $thecommand } errmsg
    puts $errmsg

    if { [ file readable $invxfmfilename ] == 0 } {
	puts "No transformation saved!"
	return 4
    }

    return 0
}

#}}}
#{{{ call GUI

wm withdraw .
invertxfm .rename
tkwait window .rename

#}}}
