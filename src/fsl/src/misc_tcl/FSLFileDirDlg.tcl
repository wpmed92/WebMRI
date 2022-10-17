# FileDirDlg.tcl --
#
#	Implements the File/Directory Selection Dialog widget.
#
# Copyright (c) 1996, Expert Interface Technologies
#
# See the file "license.terms" for information on usage and redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

tixWidgetClass FSLFileDirSelectDialog {
    -classname tixFSLFileDirSelectDialog
    -superclass tixStdDialogShell
    -method {
    }
    -flag {
	-command -filterhist -directory -pattern
    }
    -configspec {
	{-command command Command {}}
	{-title title Title "Select A File or Directory"}
	{-filterhist filterhist Filterhist {}}
	{-pattern pattern Pattern *}
	{-directory directory Directory {}}
    }
}

proc FSLFileDirSelectDialog:ConstructTopFrame {w frame} {
    upvar #0 $w data

    tixChainMethod $w ConstructTopFrame $frame

    set data(w:fsbox) [FSLFileDirSelectBox $frame.fsbox \
	-command "FSLFileDirSelectDialog:Invoke $w"\
	-filterhist $data(-filterhist) \
	-directory $data(-directory) \
	-pattern $data(-pattern)
    ]
    pack $data(w:fsbox) -expand yes -fill both
}

proc FSLFileDirSelectDialog:SetBindings {w} {
    upvar #0 $w data

    tixChainMethod $w SetBindings

    $data(w:btns) subwidget ok     config -command "$data(w:fsbox) invoke" \
	-underline 0
    $data(w:btns) subwidget apply  config -command "$data(w:fsbox) filter" \
	-text Filter -underline 0
    $data(w:btns) subwidget cancel config -command "wm withdraw $w" \
	-underline 0
    $data(w:btns) subwidget help config -underline 0


    bind $w <Alt-Key-l> "focus [$data(w:fsbox) subwidget filelist]"
    bind $w <Alt-Key-d> "focus [$data(w:fsbox) subwidget dirlist]"
    bind $w <Alt-Key-s> "focus [$data(w:fsbox) subwidget selection]"
    bind $w <Alt-Key-t> "focus [$data(w:fsbox) subwidget filter]"
    bind $w <Alt-Key-o> "tkButtonInvoke [$data(w:btns) subwidget ok]"
    bind $w <Alt-Key-f> "tkButtonInvoke [$data(w:btns) subwidget apply]"
    bind $w <Alt-Key-c> "tkButtonInvoke [$data(w:btns) subwidget cancel]"
    bind $w <Alt-Key-h> "tkButtonInvoke [$data(w:btns) subwidget help]"
}

proc FSLFileDirSelectDialog:Invoke {w filename} {

    set filename [ fix_cygwin_filename $filename ]

    upvar #0 $w data

    wm withdraw $w

    if {$data(-command) != {}} {
	set bind(specs) "%V"
	set bind(%V) $filename
	tixEvalCmdBinding $w $data(-command) bind $filename
    }
}
