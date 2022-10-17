
tixWidgetClass FSLFileEntry {
    -superclass tixPrimitive
    -classname tixFSLFileEntry
    -method {

    }
    -flag {
	-variable -directory -pattern -title -width -showdotfiles
	-dirasfile -labelwidth -command -filterhist
    }
    -configspec {
	{-variable variable Variable N}
	{-pattern pattern Pattern *}
	{-directory directory Directory {}}
	{-label label Label "File"}
	{-title title Title "FSL Select File"}
	{-width width Width 40}
	{-labelwidth labelwidth Labelwidth {}}
	{-showdotfiles showdotfiles Showdotfiles false tixVerifyBoolean}
	{-dirasfile dirasfile Dirasfile {}}
	{-command command Command {}}
	{-filterhist filterhist Filterhist {}}
    }
}

proc FSLFileEntry:ConstructWidget { w } {

    global FSLDIR

    upvar #0 $w data

    set FILEBMP ${FSLDIR}/tcl/file.xbm

    tixChainMethod $w ConstructWidget

    set data(w:fsbox) [FSLFileDirSelectDialog $w.browser \
	    -title $data(-title) \
	    -command "FSLFileEntry:Invoke $w" \
	    -filterhist $data(-filterhist) \
	    -directory $data(-directory) \
	    -pattern $data(-pattern)
    ]

    $w.browser subwidget fsbox config -showdotfiles $data(-showdotfiles)
    $w.browser subwidget fsbox config -dirasfile $data(-dirasfile)

    frame $w.frame
    set data(w:label) [ label $w.frame.label -text $data(-label) ]
    if {$data(-labelwidth) != ""} {
	$w.frame.label configure -width $data(-labelwidth)
    }

    label $w.frame.padlbl -width 1

    frame $w.frame.f -borderwidth 2 -relief sunken
    set data(w:filename) [ entry $w.frame.f.filename -textvariable $data(-variable) \
			     -width $data(-width) -borderwidth 0 -highlightthickness 0 ]

    set im [image create bitmap -file $FILEBMP]

    set data(w:button) [ button $w.frame.f.browse -text "Browse" -command "$w.browser popup" \
			     -borderwidth 2 -relief raised -pady 0 -padx 0 \
			     -image $im -highlightthickness 0 ]

    pack $w.frame.f.filename -expand yes -fill x \
	    -in $w.frame.f -side left -padx 0 -pady 0

    pack $w.frame.f.browse \
	    -in $w.frame.f -side left -padx 0 -pady 0

    pack $w.frame.label $w.frame.padlbl $w.frame.f \
	    -in $w.frame -side left

    pack $w.frame.f \
	    -in $w.frame -side left -expand yes -fill x 

    pack $w.frame -in $w -expand yes -fill x 
}

proc FSLFileEntry:Invoke { w filename } {

    set filename [ fix_cygwin_filename $filename ]

    upvar #0 $w data
    upvar $data(-variable) Variable
    
    set Variable $filename

    if {$data(-command) != {}} {
	set bind(specs) "%V"
	set bind(%V) $filename
	tixEvalCmdBinding $w $data(-command) bind $filename
    }
}

