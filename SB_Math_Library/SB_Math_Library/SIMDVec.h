/******************************************************************************/
/*!
\file SIMDVec.h
\author Javon Lee Xiong Wei, xiongweijavon.lee, 1702507
\par xiongweijavon.lee\@digipen.edu
\date 28/6/2019
\brief
This file contains the interface of the base class for SIMD Vectors

Copyright (C) 2018 DigiPen Institute of Technology.
Reproduction or disclosure of this file or its contents without the
prior written consent of DigiPen Institute of Technology is prohibited.
*******************************************************************************/
#pragma once
#include <immintrin.h> //: AVX, operations on integers, 8 float or 4 double.

//! Forward declaration for the traits class 
template <typename T>
    struct SIMD_Vec_traits;

//! This is a templated base class for all SIMD Vectors to inherit from
    template <typename T>
    class SIMD_Vec
    {
    public:
        //! A definition the value type 
        typedef typename SIMD_Vec_traits<T>::value_type value_type; 
        /*!
            \brief 
            a downcast operator to allow the user to call methods using the 
            inheriting class 
        */
        inline T& operator() ()
        {
            return static_cast<T*>(this);
        }
        //! The const version for downcasting 
        inline T& operator() () const 
        {
            return static_cast<const T*>(this);
        }

        /*! 
        The additional/subtract/multiplication/divide assignment 
        operator will call the derived classes' overloaded operators 
        */
        inline T& operator+= (const T& rhs)
        {
            (*this)() = (*this)() + rhs;
            return (*this)();
        }
        inline T& operator+= (const value_type& rhs)
        {
            (*this)() = (*this)() + static_cast<T>(rhs);
            return (*this)();
        }
        
        inline T& operator-= (const T& rhs)
        {
            (*this)() = (*this)() - rhs;
            return (*this)();
        }
        inline T& operator-= (const value_type& rhs)
        {
            (*this)() = (*this)() - static_cast<T>(rhs);
            return (*this)();
        }

        inline T& operator*= (const T& rhs)
        {
            (*this)() = (*this)() * rhs;
            return (*this)();
        }
        inline T& operator*= (const value_type& rhs)
        {
            (*this)() = (*this)() * static_cast<T>(rhs);
            return (*this)();
        }

        inline T& operator/= (const T& rhs)
        {
            (*this)() = (*this)() / rhs;
            return (*this)();
        }
        inline T& operator/= (const value_type& rhs)
        {
            (*this)() = (*this)() / static_cast<T>(rhs);
            return (*this)();
        }

        //! pre-increment operator 
        inline T& operator++()
        {
            (*this)() += value_type(1);
            return (*this)();
        }
        //! post-increment operator 
        inline T operator++(int)
        {
            T tmp = (*this)();
            (*this)() += value_type(1);
            return tmp;
        }

        //! pre-decrement operator 
        inline T& operator--()
        {
            (*this)() -= value_type(1);
            return (*this)();
        }
        //! post-decrement operator 
        inline T operator--(int)
        {
            T tmp = (*this)();
            (*this)() -= value_type(1);
            return tmp;
        }

    protected:
        /*! 
        Ensures that only the inheritting SIMD Vector classes can 
        instantiate the object 
        */
        inline  SIMD_Vec()  {}
        inline ~SIMD_Vec()  {}
        inline  SIMD_Vec(const SIMD_Vec& rhs) {}
        inline  SIMD_Vec& operator= (const SIMD_Vec& rhs) 
        {
            return *this;
        }

    };

