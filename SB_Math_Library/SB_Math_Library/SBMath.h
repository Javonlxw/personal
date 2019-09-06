/******************************************************************************/
/*!
\file SBMath.h
\author Javon Lee Xiong Wei, xiongweijavon.lee, 1702507
\par xiongweijavon.lee\@digipen.edu
\date 28/6/2019
\brief
This file contains the necessary macro definitions that will be used in various
places across SukaBread's math library, as well as the necessary includes, e.g 
the intrisincs header to allow for intrisinc operations. 

Copyright (C) 2018 DigiPen Institute of Technology.
Reproduction or disclosure of this file or its contents without the
prior written consent of DigiPen Institute of Technology is prohibited.
*******************************************************************************/
#pragma once
#include <stdint.h>    // to specify exact-width integer types
#include <math.h>      // basic math library functions 
#include "Vec3.h"  // The header for SukaBread's Vector library 

//! For commonly used functions for trigo/angle-related operations
#define PI_ 3.14159265358979323846f

//! For defining the maximum and minimum values in integer format 
#define INT_MIN     (-2147483647 - 1
#define INT_MAX     2147483647
#define FLT_MAX     3.402823466e+38F

//! For making a function inline (so that function-calling time will be saved)
#define SB_INLINE	__forceinline

//Angle Conversion
static constexpr float Deg_to_Rad(float degree) { return degree * (PI_ / 180.f); }
static constexpr float Rad_to_Deg(float radian) { return radian * (180.f / PI_); }


//For comparing values on which is bigger 
template<typename T>
T SB_MAX(T a, T b)
{
    return a > b ? a : b;
}

//For comparing values on which is smaller
template<typename T>
T SB_MIN(T a, T b)
{
    return a > b ? a : b;
}








