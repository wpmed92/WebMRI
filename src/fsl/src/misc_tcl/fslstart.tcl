
# test for whether we are in tclsh or wish
if { [ info exists env(TCLTKSHELL) ] &&  $env(TCLTKSHELL) == "wish" } {
    package require Tk
    package require Tix
    tix configure -scheme TixGray
    option add *AutoRepeat false
}

set FSLDIR $env(FSLDIR)

set FSLCONVERT $env(FSLCONVERT)
set FSLGNUPLOT $env(FSLGNUPLOT)
set FSLOUTPUTTYPE $env(FSLOUTPUTTYPE)
set FSLTCLSH $env(FSLTCLSH)
set FSLWISH $env(FSLWISH)
set FSLDISPLAY $env(FSLDISPLAY)
set FSLBROWSER $env(FSLBROWSER)

set USER       $env(USER)
set HOME       [ exec sh -c " cd ; pwd " ]
set PWD        [ exec sh -c " pwd " ]
set PROCID     [ pid ]
set HOSTNAME   [ exec hostname ]
set OS         [ exec uname ]

set auto_path [ linsert $auto_path 0 $FSLDIR/tcl ]

if { ! [ info exists INGUI ] } {
    source $FSLDIR/tcl/FSLFileDirBox.tcl
    source $FSLDIR/tcl/FSLFileDirDlg.tcl
    source $FSLDIR/tcl/FSLFileEntry.tcl
}

# what type of OS?
set OSFLAVOUR unix
set gui_ext ""
set FSLSLASH ""
set UNAME [ exec uname ]
if { $UNAME == "Darwin" } {
    set OSFLAVOUR macos
    set gui_ext "_gui"
} elseif { [ string compare CYGWIN [ string range $UNAME 0 5 ] ] == 0 } {
    set OSFLAVOUR cygwin
    set gui_ext "_gui"
    set FSLSLASH "C:/cygwin"
}

