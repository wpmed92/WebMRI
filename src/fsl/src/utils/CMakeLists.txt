PROJECT(UTILS)

FIND_PATH(FSL_DIR etc/fslconf/fsl.sh /usr/local/fsl)
INCLUDE_DIRECTORIES(${FSL_DIR}/include)
INCLUDE_DIRECTORIES(${FSL_DIR}/extras/include/newmat)

IF(NOT LIBRARY_OUTPUT_PATH)
  SET(LIBRARY_OUTPUT_PATH ${UTILS_BINARY_DIR}/lib)
ENDIF(NOT LIBRARY_OUTPUT_PATH)

SET(libsrcs
  matches.cc functions.cc usage.cc check.cc parse.cc log.cc time_tracer.cc)

ADD_LIBRARY(utils ${libsrcs})
ADD_EXECUTABLE(opttst opttst.cc)
ADD_DEPENDENCIES(opttst utils)

TARGET_LINK_LIBRARIES(opttst ${LIBRARY_OUTPUT_PATH}/libutils.a)

INSTALL_FILES(/include/utils .*\\.h$ )
INSTALL_TARGETS(/lib utils)

ENABLE_TESTING()
ADD_TEST(ShouldPass opttst -d 7 -m mask -r res)
ADD_TEST(ShouldFail opttst -V)
ADD_TEST(ShouldFail opttst -verb)
