# FSL configuration file 
#  - to be sourced by the user, typically in .bashrc or equivalent
#  - note that the user should set 

# Written by Mark Jenkinson
#  FMRIB Analysis Group, University of Oxford

# SHCOPYRIGHT


#### Set up standard FSL user environment variables ####

# The following variable selects the default output image type
# Legal values are:  ANALYZE  NIFTI  NIFTI_PAIR  ANALYZE_GZ  NIFTI_GZ  NIFTI_PAIR_GZ
# This would typically be overwritten in ${HOME}/.fslconf/fsl.sh if the user wished
#  to write files with a different format
FSLOUTPUTTYPE=NIFTI
export FSLOUTPUTTYPE

# Comment out the definition of FSLMULTIFILEQUIT to enable 
#  FSL programs to soldier on after detecting multiple image
#  files with the same basename ( e.g. epi.hdr and epi.nii )
FSLMULTIFILEQUIT=TRUE ; export FSLMULTIFILEQUIT


# The following variables specify paths for programs and can be changed
#  or replaced by different programs ( e.g. FSLDISPLAY=open   for MacOSX)

FSLTCLSH=$FSLDIR/bin/tclsh
FSLWISH=$FSLDIR/bin/wish
FSLGNUPLOT=$FSLDIR/bin/gnuplot
FSLDISPLAY=$FSLDIR/bin/display
FSLCONVERT=$FSLDIR/bin/convert

FSLBROWSER=$FSLDIR/tcl/fslwebbrowser

export FSLTCLSH FSLWISH FSLGNUPLOT FSLDISPLAY FSLCONVERT FSLBROWSER


# The following variables are used for running code in parallel across
#  several machines ( i.e. for FDT )

FSLLOCKDIR=
FSLMACHINELIST=
FSLREMOTECALL=

export FSLLOCKDIR FSLMACHINELIST FSLREMOTECALL

# Set up development variables (not for the faint-hearted)

FSLCONFDIR=$FSLDIR/config
FSLMACHTYPE=
if [ "`which gcc | grep no`" = "" ] ; then
    FSLMACHTYPE=`gcc -dumpmachine`-gcc`gcc -dumpversion`;
    if [ `uname` = Darwin ] ; then
	if [ "`system_profiler -detailLevel -2 | grep -i 'Machine Name' | grep G5`" != "" ] ; then
	    FSLMACHTYPE=${FSLMACHTYPE}-G5
	fi
    fi
fi
export FSLCONFDIR FSLMACHTYPE


###################################################
### Add other global environment variables here ###
###      or change the definitions above        ###
###################################################


# USER VARIABLES HERE


###################################################
####    DO NOT ADD ANYTHING BELOW THIS LINE    ####
###################################################

if [ -f /usr/local/etc/fslconf/fsl.sh ] ; then
  . /usr/local/etc/fslconf/fsl.sh ;
fi


if [ -f /etc/fslconf/fsl.sh ] ; then
  . /etc/fslconf/fsl.sh ;
fi


if [ -f "${HOME}/.fslconf/fsl.sh" ] ; then
  . "${HOME}/.fslconf/fsl.sh" ;
fi
