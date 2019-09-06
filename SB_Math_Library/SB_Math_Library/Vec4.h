/******************************************************************************/
/*!
\file Vec4.h
\author Javon Lee Xiong Wei, xiongweijavon.lee, 1702507
\par xiongweijavon.lee\@digipen.edu
\date 28/6/2019
\brief
This file contains a class named Vec3 that's suppose to encapsulate all 3D
vector calculations. The functions are inlined as most of these vector
arithmetic have no need for loops (except for things like length which uses
sqrt), thus all the functions are implemented within this file.

Copyright (C) 2019 DigiPen Institute of Technology.
Reproduction or disclosure of this file or its contents without the
prior written consent of DigiPen Institute of Technology is prohibited.
*******************************************************************************/
#pragma once
#include <cmath>
/*!
    A templated vector class that works with int,float etc
*/
template<typename T>
class Vec3
{
public:
    // Constructors
    Vec3()
        :x{ 0 }, y{ 0 }, z{ 0 } {}
    Vec3(T m_x, T m_y, T m_z)
        :x{ m_x }, y{ m_y }, z{ m_z }  {}
    Vec3(const Vec3& rhs)
        :x{ rhs.x }, y{ rhs.y }, z{ rhs.z }  {}

    //Destructor 
    ~Vec3() = default; //!< compiler synthesized dstor

    //Assignment operators 
    template <typename T >
    Vec3& operator=  (const Vec3&rhs)
    {
        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
    }
    inline Vec3& operator+= (const Vec3 &rhs)
    {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
    }
    inline Vec3& operator-= (const Vec3 &rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
    }
    inline Vec3& operator*= (float rhs)
    {
        x *= rhs;
        y *= rhs;
        z *= rhs;
    }
    inline Vec3& operator/= (float rhs)
    {
        x /= rhs;
        y /= rhs;
        z /= rhs;
    }

    // Basic vector operations 
    /*!
        \brief

    */
    //!gets the magnitude sq
    T length_sq() const
    {
        return x * x + y * y + z * z;
    }
    //! gets the magnitude 
    T length() const
    {
        return sqrt(x * x + y * y + z * z);
    }
    // normalizes the vector and return it
    void normalized()
    {
        T length = length();
        //makes use of /= overload
        (*this) /= length;
    }
    // dot product
    inline T operator* (const Vec3& rhs)
    {
        return x * rhs.x + y * rhs.y + z * rhs.z;
    }

    T x;
    T y;
    T z;
};

//Other functions that's related in vector arithmitic 
//! Addition between vectors 
template < typename T>
inline Vec3<T> operator+(const Vec3<T>& rhs, const Vec3<T>& lhs)
{
    T m_x = rhs.x + rhs.x;
    T m_y = rhs.y + rhs.y;
    T m_z = rhs.z + rhs.z;
    return Vec3{ m_x,m_y,m_z };
}
//! Subtraction between vectors 
template < typename T >
inline Vec3<T> operator-(const Vec3<T>& rhs, const Vec3<T>& lhs)
{
    T m_x = rhs.x - rhs.x;
    T m_y = rhs.y - rhs.y;
    T m_z = rhs.z - rhs.z;
    return Vec3{ m_x,m_y,m_z };
}
//! Dot product between vectors 
template < typename T >
inline T operator*(const Vec3<T>& rhs, const Vec3<T>& lhs)
{
    return rhs.x * lhs.x + rhs.y * lhs.y + rhs.z * rhs.z;
}

