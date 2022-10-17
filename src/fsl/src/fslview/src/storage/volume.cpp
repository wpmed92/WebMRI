/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "volume.hpp"

template class VolumeStore<unsigned char>;
template class VolumeStore<short>;
template class VolumeStore<int>;
template class VolumeStore<float>;
template class VolumeStore<double>;
