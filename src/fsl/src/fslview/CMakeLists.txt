PROJECT(FSLVIEW)

SET(LIBRARY_OUTPUT_PATH ${FSLVIEW_BINARY_DIR}/lib)
SET(EXECUTABLE_OUTPUT_PATH ${FSLVIEW_BINARY_DIR}/bin)

FIND_PATH( FSL_DIR etc/fslconf/fsl.sh
           /usr/local/fsl )
INCLUDE_DIRECTORIES(${FSL_DIR}/src)
INCLUDE_DIRECTORIES(${FSL_DIR}/extras/src/newmat)
INCLUDE_DIRECTORIES( ${FSL_DIR}/extras/include/zlib )

INCLUDE_DIRECTORIES(${FSLVIEW_SOURCE_DIR}/src)

SUBDIRS(${FSL_DIR}/src/fslio ${FSL_DIR}/src/niftiio ${FSL_DIR}/src/znzlib)
SUBDIRS(${FSL_DIR}/src/miscmaths ${FSL_DIR}/extras/src/newmat)
SUBDIRS(src)
