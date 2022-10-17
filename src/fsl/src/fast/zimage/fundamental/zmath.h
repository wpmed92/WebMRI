/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file defines some basic mathematical constants and should be included
// as the single headfile for matrix/vector calculations.

#ifndef __INC_ZMATH_H__
#define __INC_ZMATH_H__
#include <cmath>

#undef M_E
#define M_E				2.7182818284590452354	// e
#undef M_LOG2E
#define M_LOG2E			1.4426950408889634074	// log2(e)
#undef M_LOG10E
#define M_LOG10E		0.43429448190325182765	// log10(e)
#undef M_LN2
#define M_LN2			0.69314718055994530942	// ln(2)
#undef M_LN10
#define M_LN10			2.3025850929940456840	// ln(10)

#undef M_PI
#define M_PI			3.1415926535897932385	// pi
#undef M_2PI
#define M_2PI			6.283185307179586477	// 2 * pi
#undef M_PI_2
#define M_PI_2			1.5707963267948966192	// pi/2
#undef M_PI_4
#define M_PI_4			0.78539816339744830962	// pi/4
#undef M_1_PI
#define M_1_PI			0.31830988618379067154	// 1/pi
#undef M_2_PI
#define M_2_PI			0.636619772367581343	// 2/pi
#undef M_SQRT_PI
#define M_SQRT_PI		1.772453850905516		// sqrt(pi)
#undef M_2_SQRTPI
#define M_2_SQRTPI		1.1283791670955125739	// 2/sqrt(pi)
#undef M_1_SQRT_2PI
#define M_1_SQRT_2PI	0.398942280401432678	// 1/sqrt(2*PI)
#undef M_SQRT2
#define M_SQRT2			1.4142135623730950488	// sqrt(2)
#undef M_SQRT1_2
#define M_SQRT1_2		0.7071067811865475244	// 1/sqrt(2)
#undef M_SQRT3
#define M_SQRT3			1.73205080756887729		// sqrt(3)
#undef M_SQRT1_3
#define M_SQRT1_3		0.6933612743506347		// 1/sqrt(3)

#include "mathcommon.h"
#include "eigensystem.h"
#include "svd.h"

#endif // __INC_ZMATH_H__
