/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "operator.h"
#include "mydefine.h"

#ifdef _MSC_VER

template <> double Opt<char>::min = -128; //MinValue<char>();
template <> double Opt<char>::max = 127; //MaxValue<char>();
template <> double Opt<BYTE>::min = 0; //MinValue<BYTE>();
template <> double Opt<BYTE>::max = 255; //MaxValue<BYTE>();
template <> double Opt<short>::min = -32768.0; //MinValue<short>();
template <> double Opt<short>::max = 32767; //MaxValue<short>();
template <> double Opt<WORD>::min = 0; //MinValue<WORD>();
template <> double Opt<WORD>::max = 65536; // MaxValue<WORD>();
template <> double Opt<int>::min = -2147483648.0; //MinValue<int>();
template <> double Opt<int>::max = 2147483647.0; //MaxValue<int>();
template <> double Opt<UINT>::min = 0; //MinValue<UINT>();
template <> double Opt<UINT>::max = 4294967295.0; //MaxValue<UINT>();
template <> double Opt<float>::min = -2147483648.0; //MinValue<float>();
template <> double Opt<float>::max = 2147483647.0; //MaxValue<float>();
template <> double Opt<double>::min = -2147483648.0; //MinValue<double>();
template <> double Opt<double>::max = 2147483647.0; //MaxValue<double>();

#endif
