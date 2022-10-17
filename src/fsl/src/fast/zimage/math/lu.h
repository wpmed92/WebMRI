/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file declares functions for LU decomposing.

#ifndef __LU_H_  // prevent multiple includes
#define __LU_H_

bool LUdcmp(ZMatrix<double>& A);
bool LUdcmp(ZMatrix<double>& A, ZVector<int>& indx);
bool LUdcmp(ZMatrix<float>& A);
bool LUdcmp(ZMatrix<float>& A, ZVector<int>& indx);
void LUbksb (const ZMatrix<double>& A, const ZVector<int>& indx, ZVector<double>& B );
void LUbksb (const ZMatrix<float>& A, const ZVector<int>& indx, ZVector<float>& B );

#endif
