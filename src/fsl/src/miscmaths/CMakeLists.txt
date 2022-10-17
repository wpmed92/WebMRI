SET(sources
  miscmaths.cc
  optimise.cc
  miscprob.cc
  kernel.cc
  histogram.cc
  base2z.cc
  t2z.cc
  f2z.cc
  volume.cc
  volumeseries.cc
  minimize.cc
  cspline.cc
  sparse_matrix.cc
  sparsefn.cc
  rungekutta.cc)

ADD_LIBRARY(miscmaths ${sources})

INCLUDE_DIRECTORIES(${FSL_DIR}/include)
INCLUDE_DIRECTORIES(${FSL_DIR}/extras/include/newmat)
INCLUDE_DIRECTORIES(${FSL_DIR}/extras/include/libprob)

INSTALL_FILES(/include/miscmaths .h ${sources})