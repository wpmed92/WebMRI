# FileDirBox.tcl --
#
#	Implements the Folder Selection Box widget.
#
# Copyright (c) 1996, Expert Interface Technologies
#
# See the file "license.terms" for information on usage and redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.
#


# ToDo
#   (1)	If user has entered an invalid directory, give an error dialog
#

tixWidgetClass FSLFileDirSelectBox {
    -superclass tixPrimitive
    -classname  tixFSLFileDirSelectBox
    -method {
	filter invoke
    }
    -flag {
	-browsecmd -command -dir -directory -disablecallback
	-grab -pattern -selection -value -showdotfiles
	-dirasfile -filterhist
    }
    -configspec {
	{-multiplefiles multipleFiles MultipleFiles false tixVerifyBoolean }
	{-browsecmd browseCmd BrowseCmd {}}
	{-command command Command {}}
	{-directory directory Directory {}}
	{-disablecallback disableCallback DisableCallback false tixVerifyBoolean}
	{-showdotfiles showdotfiles Showdotfiles false tixVerifyBoolean}
	{-grab grab Grab global}
	{-pattern pattern Pattern *}
	{-value value Value {}}
	{-dirasfile dirasfile Dirasfile {}}
	{-filterhist filterhist Filterhist {}}
    }
    -alias {
	{-selection -value}
	{-dir -directory}
    }
    -forcecall {
	-value
    }
    -default {
	{.relief			raised}
	{*filelist*Listbox.takeFocus	true}
	{.borderWidth 			1}
	{*Label.anchor			w}
	{*Label.borderWidth		0}
	{*Label.font                   -Adobe-Helvetica-Bold-R-Normal--*-120-*}
	{*TixComboBox*scrollbar		auto}
	{*TixComboBox*Label.anchor	w}
	{*TixScrolledListBox.scrollbar	auto}
	{*Listbox.exportSelection	false}
	{*Listbox.exportSelection	false}
    }
}

proc FSLFileDirSelectBox:InitWidgetRec {w} {
    upvar #0 $w data
    global env

    tixChainMethod $w InitWidgetRec

    if {$data(-directory) == {}} {
	global env

	if {[info exists env(PWD)]} {
	    set data(-directory) $env(PWD)
	} else {
	    set data(-directory) [pwd]
	}
    }

    set data(flag)      0
    set data(fakeDir)   0
}

#----------------------------------------------------------------------
#		Construct widget
#----------------------------------------------------------------------
proc FSLFileDirSelectBox:ConstructWidget {w} {
    upvar #0 $w data

    tixChainMethod $w ConstructWidget

    set frame1 [FSLFileDirSelectBox:CreateFrame1 $w]
    set frame2 [FSLFileDirSelectBox:CreateFrame2 $w]
    set frame3 [FSLFileDirSelectBox:CreateFrame3 $w]

    pack $frame1 -in $w -side top -fill x
    pack $frame3 -in $w -side bottom -fill x
    pack $frame2 -in $w -side top -fill both -expand yes

    tixSetSilent $data(w:filter) \
	[FSLFileDirSelectBox:GetFilter $w $data(-directory) $data(-pattern)]

    set filter [FSLFileDirSelectBox:GetFilter $w $data(-directory) $data(-pattern)]
    
    $data(w:filter) addhistory \
	    [FSLFileDirSelectBox:InterpFilter $w $filter]
}

proc FSLFileDirSelectBox:CreateFrame1 {w} {
    upvar #0 $w data

    frame $w.f1 -border 10
    tixComboBox $w.f1.filter -history true\
	-command "$w filter" -anchor e \
	-label "Filter:" \
	-options {
	    slistbox.scrollbar auto
	    listbox.height 5
	    label.anchor w
            label.underline 3
            label.side top
	}
    # External filter history
    if {$data(-filterhist) != {}} {
	$w.f1.filter configure \
		-listcmd "FSLFileDirSelectBox:FilterHist $w"
    }

    set data(w:filter) $w.f1.filter

    pack $data(w:filter) -side top -expand yes -fill both
    return $w.f1
}

proc FSLFileDirSelectBox:CreateFrame2 {w} {
    upvar #0 $w data

    tixPanedWindow $w.f2 -orientation horizontal -height 200

    #     THE LEFT FRAME
    #-----------------------
    set dir [$w.f2 add directory -size 200]
    $dir config -relief flat
    label $dir.lab -text "Directories:" -underline 0
    set data(w:dirlist) [tixScrolledListBox $dir.dirlist \
		       -scrollbar auto \
		       -options { listbox.width 4 listbox.height 6 }]

    $data(w:dirlist) subwidget listbox config  -selectmode single


    pack $dir.lab -side top -fill x -padx 10
    pack $data(w:dirlist) -side bottom -expand yes -fill both -padx 10

    #     THE RIGHT FRAME
    #-----------------------
    set file [$w.f2 add file -size 200]
    $file config -relief flat
    label $file.lab -text "Files:" -underline 2
    set data(w:filelist) [tixScrolledListBox $file.filelist \
		       -scrollbar auto\
		       -options {listbox.width 4 listbox.height 6}]

    pack $file.lab -side top -fill x -padx 10
    pack $data(w:filelist) -side bottom -expand yes -fill both -padx 10

    return $w.f2
}

proc FSLFileDirSelectBox:CreateFrame3 {w} {
    upvar #0 $w data

    frame $w.f3 -border 10
    tixComboBox $w.f3.selection -history true\
	-command "FSLFileDirSelectBox:SelInvoke $w" \
	-anchor e \
	-label "Selection:" \
	-options {
	    slistbox.scrollbar auto
	    listbox.height 5
	    label.anchor w
            label.underline 0
            label.side top
	}

    set data(w:selection) $w.f3.selection

    pack $data(w:selection) -side top -fill both

    return $w.f3
}

proc FSLFileDirSelectBox:SelInvoke {w args} {
    upvar #0 $w data

    set event [tixEvent type]

    $data(w:filelist) subwidget listbox selection clear 0 end

    if {$event != "<FocusOut>" && $event != "<Tab>"} {
	$w invoke
    }
}


#----------------------------------------------------------------------
#                           BINDINGS
#----------------------------------------------------------------------

proc FSLFileDirSelectBox:SetBindings {w} {
    upvar #0 $w data

    tixChainMethod $w SetBindings

    tixDoWhenMapped $w "FSLFileDirSelectBox:FirstMapped $w"

    $data(w:dirlist) config \
	-browsecmd "FSLFileDirSelectBox:SelectDir $w" \
	-command   "FSLFileDirSelectBox:InvokeDir $w"

    $data(w:filelist) config \
	-browsecmd "FSLFileDirSelectBox:SelectFile $w" \
	-command   "FSLFileDirSelectBox:InvokeFileDir $w"

    set listbox [ $data(w:filelist) subwidget listbox ]
    if { $data(-multiplefiles) } {
	$listbox config -selectmode extended
    } else {
	$listbox config -selectmode single
    }

}

#----------------------------------------------------------------------
#                           CONFIG OPTIONS
#----------------------------------------------------------------------
proc FSLFileDirSelectBox:config-directory {w value} {
    upvar #0 $w data

    set value [tixFile tildesubst $value]
    set value [tixFile trimslash $value]

    tixSetSilent $data(w:filter) \
	[FSLFileDirSelectBox:GetFilter $w $value $data(-pattern)]

    $w filter
    set data(-directory) $value
    return $value
}

proc FSLFileDirSelectBox:config-pattern {w value} {
    upvar #0 $w data

    if {$value == {}} {
	set data(-pattern) "*"
    } else {
	set data(-pattern) $value
    }
    
    tixSetSilent $data(w:filter) \
	[FSLFileDirSelectBox:GetFilter $w $data(-directory) $value]

    # Returning a value means we have overridden the value and updated
    # the widget record ourselves.
    #
    return $data(-pattern)
}

proc FSLFileDirSelectBox:config-value {w value} {
    upvar #0 $w data

    tixSetSilent $data(w:selection) $value
}

#----------------------------------------------------------------------
#                    PUBLIC METHODS
#----------------------------------------------------------------------
proc FSLFileDirSelectBox:filter {w args} {
    upvar #0 $w data

    $data(w:filter) popdown
    set filter [FSLFileDirSelectBox:InterpFilter $w]
    FSLFileDirSelectBox:LoadDir $w
}

# InterpFilter:
#	Interp the value of the w:filter widget. 
#
# Side effects:
#	Changes the fields data(-directory) and data(-pattern) 
#
proc FSLFileDirSelectBox:InterpFilter {w {filter {}}} {
    upvar #0 $w data

    if {$filter == {}} {
	set filter [$data(w:filter) cget -selection]
	if {$filter == {}} {
	    set filter [$data(w:filter) cget -value]
	}
    }

    set filter [tixFile tildesubst $filter]
    
    if [file isdir $filter] {
	set data(-directory) [tixFile trimslash $filter]
	set data(-pattern) "*"
    } else {
	set data(-directory)  [file dir $filter]
	set data(-pattern)    [file tail $filter]
    }

    set OWD [ pwd ]
    cd $data(-directory)
    set data(-directory) [ pwd ]
    cd $OWD

    set filter [FSLFileDirSelectBox:GetFilter $w $data(-directory) \
        $data(-pattern)]

    tixSetSilent $data(w:filter) $filter

    return $filter
}

proc FSLFileDirSelectBox:invoke {w args} {
    upvar #0 $w data
    upvar #0 $data(-filterhist) FilterHist

    if {[$data(w:selection) cget -value] !=
	[$data(w:selection) cget -selection]} {
	    $data(w:selection) invoke
	    return
    }
    
    # record the filter
    #
    set filter [FSLFileDirSelectBox:InterpFilter $w]
    $data(w:filter) addhistory $filter

    # External filter history
    if {$data(-filterhist) != {}} {
	lappend FilterHist [file dir $filter]
    }

    # record the selection
    #
    set value [$data(w:selection) cget -value]
    set value [tixFile tildesubst $value]
    set value [tixFile trimslash $value]

    if {[string index $value 0] != "/"} {
	set value $data(-directory)/$value
	set value [tixFile tildesubst $value]
	set value [tixFile trimslash $value]
	tixSetSilent $data(w:selection) $value
    }
    set data(-value) $value

    $data(w:selection) addhistory $data(-value)

    $data(w:filter) align
    $data(w:selection) align

    if {$data(-command) != {} && !$data(-disablecallback)} {
	set bind(specs) "%V"

	if { $data(-multiplefiles) } {
	    set listbox  [ $data(w:filelist) subwidget listbox ] 
	    set lfiles [ $listbox curselection ]


	    if { [ llength $lfiles ] == 0 } {
		set sfile [ list $data(-value) ]
	    } else {
		set sfile ""
		foreach i $lfiles {
		    set i [ $listbox get $i ]

		    set i [tixFile tildesubst $i]
		    set i [tixFile trimslash $i]

		    if {[string index $i 0] != "/"} {
			set i $data(-directory)/$i
			set i [tixFile tildesubst $i]
			set i [tixFile trimslash $i]
		    }

		    lappend sfile $i
		}
		if { [ llength $sfile ] == 1} {
		    $data(w:selection) addhistory [ lindex $sfile 0 ]
		}
	    }


	} else {
	    set sfile $data(-value)
	}
	set bind(%V) $sfile 
	tixEvalCmdBinding $w $data(-command) bind $sfile
    }
}

#----------------------------------------------------------------------
#                    INTERNAL METHODS
#----------------------------------------------------------------------

#
# The reason for this stupid thing is that the various tcl functions
# blow up if you do "//" as they think this is some sort of remote
# file system escape.
#
proc FSLFileDirSelectBox:GetFilter {w dir pattern} {
    if {$dir == "/"} {
	return /$pattern
    } else {
	return $dir/$pattern
    }
}

proc FSLFileDirSelectBox:LoadDirIntoLists {w} {
    upvar #0 $w data

    $data(w:dirlist) subwidget listbox delete 0 end
    $data(w:filelist) subwidget listbox delete 0 end

    set appPWD [pwd]

    if [catch {cd $data(-directory)} err] {
	# The user has entered an invalid directory
	# %% todo: prompt error, go back to last succeed directory
	cd $appPWD
	return
    }

    if { $data(-directory) == "/" } {
	set dir ""
    } else {
	set dir $data(-directory)
    }

    if { $data(-showdotfiles) } {
	set dirpatmatch ".*"
    } else {
	set dirpatmatch ".."
    }

    foreach fname [lsort [glob -nocomplain * $dirpatmatch ]] {
	if {![string compare . $fname]} {
	    continue
	}

	if [file isdirectory $dir/$fname] {

	    if {$data(-dirasfile) != ""} {
		if [file isfile $dir/$fname/$data(-dirasfile)] {
		    $data(w:filelist) subwidget listbox insert end $fname
		    $data(w:dirlist) subwidget listbox insert end $fname
		} else {
		    $data(w:dirlist) subwidget listbox insert end $fname
		}
	    } else {
		$data(w:dirlist) subwidget listbox insert end $fname
	    }
	}
    }

    if {$data(-dirasfile) == ""} {
	
	if { $data(-showdotfiles) } {
	    set filepatmatch ".*"
	} else {
	    set filepatmatch ".."
	}
	
	set top 0
	if {$data(-pattern) == "*"} {
	    foreach fname [lsort [glob -nocomplain * $filepatmatch]] {
		if {![file isdirectory $dir/$fname]} {
		    $data(w:filelist) subwidget listbox insert end $fname
		    if [string match .* $fname] {
			incr top
		    }
		}
	    }
	} elseif {$data(-pattern) == "IMAGE"} {
	    foreach fname [lsort [glob -nocomplain *.hdr* *.nii* *.mnc* ]] {
		if {![file isdirectory $dir/$fname]} {
		    $data(w:filelist) subwidget listbox insert end $fname
		    if [string match .* $fname] {
			incr top
		    }
		}
	    }
	} else {
	    foreach fname [lsort [glob -nocomplain $data(-pattern)]] {
		if {![file isdirectory $dir/$fname]} {
		    $data(w:filelist) subwidget listbox insert end $fname
		}
	    }
	}
    }

    $data(w:filelist) subwidget listbox yview $top
    cd $appPWD
}

proc FSLFileDirSelectBox:LoadDir {w} {
    upvar #0 $w data

    tixBusy $w on [$data(w:dirlist) subwidget listbox]

    catch {
	# This will fail if a directory is not readable .... or some
	# strange reasons
	#
	FSLFileDirSelectBox:LoadDirIntoLists $w
	FSLFileDirSelectBox:MkDirMenu $w
    } err

    if {[$data(w:dirlist) subwidget listbox size] == 0} {
	$data(w:dirlist) subwidget listbox insert 0 ".."
    }


    tixWidgetDoWhenIdle tixBusy $w off [$data(w:dirlist) subwidget listbox]

#    if {$err != {}} {
#	error $err
#    }
}

# %% unimplemented
#
proc FSLFileDirSelectBox:MkDirMenu {w} {
    upvar #0 $w data
}

proc FSLFileDirSelectBox:SelectDir {w} {
    upvar #0 $w data

    if {$data(fakeDir) > 0} {
	incr data(fakeDir) -1
	$data(w:dirlist) subwidget listbox select clear 0 end
	$data(w:dirlist) subwidget listbox activate -1
	return
    }

    if {$data(flag)} {
	return
    }
    set data(flag) 1

    set subdir [tixListboxGetCurrent [$data(w:dirlist) subwidget listbox]]
    if {$subdir == {}} {
	set subdir "."
    }

    set filter \
	[FSLFileDirSelectBox:GetFilter $w $data(-directory) \
         $subdir/$data(-pattern)]

    tixSetSilent $data(w:filter) $filter
    
    set data(flag) 0

    # Now select the directory
    #
    set selected [tixListboxGetCurrent [$data(w:dirlist) subwidget listbox]]
    if {$selected  != {}} {
	# Make sure that the selection is not empty!
	#
	if {$data(-directory) == "/"} {
	    tixSetSilent $data(w:selection) /$selected
	    set data(-value) /$selected
	} else {
	    tixSetSilent $data(w:selection) $data(-directory)/$selected
	    set data(-value) $data(-directory)/$selected
	}
	if {$data(-browsecmd) != {}} {
	    tixEvalCmdBinding $w $data(-browsecmd) {} \
		[$data(w:selection) cget -value]
	}
    }
    set data(flag) 0
}

proc FSLFileDirSelectBox:InvokeDir {w} {
    upvar #0 $w data

    FSLFileDirSelectBox:SelectDir $w

    set theDir [$data(w:dirlist) subwidget listbox get active]

    set OWD [ pwd ]
    cd $data(-directory)/$theDir
    set data(-directory) [ pwd ]
    cd $OWD

    $data(w:dirlist) subwidget listbox select clear 0 end

    FSLFileDirSelectBox:InterpFilter $w \
	[FSLFileDirSelectBox:GetFilter $w $data(-directory) $data(-pattern)]

    FSLFileDirSelectBox:LoadDir $w

    if {![tixEvent match <Return>]} {
	incr data(fakeDir) 1
    }
}

proc FSLFileDirSelectBox:SelectFile {w} {
    upvar #0 $w data

    if {$data(flag)} {
	return
    }
    set data(flag) 1

    # Reset the "Filter:" box to the current directory:
    #	
    $data(w:dirlist) subwidget listbox select clear 0 end
    set filter \
	[FSLFileDirSelectBox:GetFilter $w $data(-directory) \
         $data(-pattern)]

    tixSetSilent $data(w:filter) $filter

    # Now select the file
    #
    set selected [tixListboxGetCurrent [$data(w:filelist) subwidget listbox]]
    if {$selected  != {}} {
	# Make sure that the selection is not empty!
	#
	if {$data(-directory) == "/"} {
	    tixSetSilent $data(w:selection) /$selected
	    set data(-value) /$selected
	} else {
	    tixSetSilent $data(w:selection) $data(-directory)/$selected
	    set data(-value) $data(-directory)/$selected
	}
	if {$data(-browsecmd) != {}} {
	    tixEvalCmdBinding $w $data(-browsecmd) {} \
		[$data(w:selection) cget -value]
	}
    }
    set data(flag) 0
}

proc FSLFileDirSelectBox:InvokeFileDir {w} {
    upvar #0 $w data

    set selected [tixListboxGetCurrent [$data(w:filelist) subwidget listbox]]
    if {$selected  != {}} {
	$w invoke
    }
}

# This is only called the first this FileDirBox is mapped -- load the directory
#
proc FSLFileDirSelectBox:FirstMapped {w} {
    if {![winfo exists $w]} {
	return
    }

    upvar #0 $w data

    FSLFileDirSelectBox:LoadDir $w
    $data(w:filter) align
}

# External filter history
proc FSLFileDirSelectBox:FilterHist { w } {

    upvar #0 $w data
    upvar #0 $data(-filterhist) FilterHist

    foreach item $FilterHist {
	set filter [FSLFileDirSelectBox:GetFilter $w $item $data(-pattern)]
	$data(w:filter) addhistory $filter
    }
}
