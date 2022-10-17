#{{{ copyright and setup

#   featwatcher - GUI for watching FEAT analysis progress
#
#   Stephen Smith, FMRIB Image Analysis Group
#
#   Copyright (C) 2003 University of Oxford
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

#}}}
#{{{ featwatcher

proc featwatcher { w } {
    global fmri FSLDIR OSFLAVOUR argc argv FSLSLASH

    toplevel $w
    wm title      $w "FEAT Watcher"
    wm iconname   $w "FEAT Watcher"
    wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm
#   wm geometry $w 800x500

    #{{{ process argument

if { [ lindex $argv 0 ] == "" } {
    puts "Usage: Featwatcher <feat_directory.feat>"
    exit 1
}

set FD [ lindex $argv 0 ]

if { ( $OSFLAVOUR != "cygwin" && [ string first / $FD ] != 0 ) ||
     ( $OSFLAVOUR == "cygwin" && [ string first : $FD ] != 1 ) } {
    set FD [ pwd ]/$FD
}

cd $FD

#}}}
    #{{{ FEAT directory name

frame $w.featdir -relief raised -borderwidth 1

label $w.featdir.label -text "$FD"

pack $w.featdir.label -in $w.featdir -padx 5 -pady 5 -side left -expand yes

#}}}
    #{{{ current command

frame $w.current -relief raised -borderwidth 1

message $w.current.label -text "" -width 780

pack $w.current.label -in $w.current -padx 5 -pady 5 -side left -expand yes

#}}}
    #{{{ log

frame $w.log -relief raised -borderwidth 1

scrollbar $w.log.sbar -command "$w.log.text yview"
text $w.log.text -yscrollcommand "$w.log.sbar set"

pack $w.log.sbar -side right -fill y
pack $w.log.text -side left -expand yes -fill both

#}}}
    #{{{ bottom buttons

frame $w.btns -relief raised -borderwidth 1

set helpfile ${FD}/report.html
if { [ string first $FSLSLASH $helpfile ] == 0 } {
    set helpfile [ string range $helpfile 9 end ]
}
button $w.btns.reportview   -text "View FEAT report webpage" -command "FmribWebHelp file: $helpfile"

button $w.btns.cancel -text "Exit Featwatcher (FEAT will continue)" -command "exit"

button $w.btns.help   -text "Help" -command "FmribWebHelp file: ${FSLDIR}/doc/feat5/featwatcher.html"

pack $w.btns.reportview $w.btns.cancel $w.btns.help -in $w.btns -padx 5 -pady 5 -side left -expand yes

#}}}

    pack $w.featdir $w.current $w.log -in $w -padx 5 -pady 5 -fill x 
    pack $w.btns -in $w -side bottom -padx 5 -pady 5 -fill x 

    set i 40
    while { 1 } {

	if { $i == 80 } {
	    if { ( $OSFLAVOUR == "unix" &&
		   [ catch { exec sh -c "ps o comm,pcpu,cputime,vsz | grep `cat .fslcurrent` | awk '{print \$1 \"   %cpu=\"\$2 \"   cputime=\"\$3 \"   mem=\"\$4}' " } fmricurrent ] == 0 ) ||
		 ( $OSFLAVOUR == "macos" &&
		   [ catch { exec sh -c "ps co command,pcpu,cputime,vsz | grep `cat .fslcurrent` | awk '{print \$1 \"   %cpu=\"\$2 \"   cputime=\"\$3 \"   mem=\"\$4}' " } fmricurrent ] == 0 ) ||
		 ( $OSFLAVOUR == "cygwin" &&
		   [ catch { exec sh -c "top n 1 b | grep `cat .fslcurrent` | awk '{print \$12 \"   %cpu=\"\$9 \"   cputime=\"\$11 \"   mem=\"\$6}' " } fmricurrent ] == 0 ) } {
		$w.current.label configure -text $fmricurrent
	    } else {
		$w.current.label configure -text ""
	    }

	    catch { exec sh -c "cat report.log" } fmrilog
	    $w.log.text delete 0.0 end
	    $w.log.text insert end $fmrilog
	    $w.log.text see end

	    set i 0
	}
	
	update
	incr i 1
	after 250
    }
}

#}}}
#{{{ call GUI and wait

if { ! [ info exists INGUI ] } {
    wm withdraw .
    featwatcher .r
    tkwait window .r
}

#}}}
