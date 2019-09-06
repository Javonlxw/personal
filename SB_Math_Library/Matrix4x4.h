#ifndef Matrix4X4_H_
#define Matrix4X4_H_

#pragma once
#include "Vector2D.h"
#include "glm\glm.hpp"
#include "glm\gtc\matrix_transform.hpp"
/**************************************************************************/
/*!

*/
/**************************************************************************/
typedef union Matrix4x4
{
  float m[4][4];

  Matrix4x4() {}
  Matrix4x4(const float *pArr);
  Matrix4x4(float _00, float _01, float _02, float _03,
    float _10, float _11, float _12, float _13,
    float _20, float _21, float _22, float _23,
    float _30, float _31, float _32, float _33);

  Matrix4x4(const glm::mat4&);
  Matrix4x4(const Matrix4x4&);
  Matrix4x4(float);

  glm::mat4 ConvertMtx44toGlmMat4();

  Matrix4x4& operator=(const Matrix4x4 &rhs);

  // Assignment operators
  Matrix4x4& operator *= (const Matrix4x4 &rhs);

} Matrix4x4, Mtx44;

Matrix4x4 operator * (const Matrix4x4 &lhs, const Matrix4x4 &rhs);

/**************************************************************************/
/*!
This operator multiplies the matrix pMtx with the vector rhs
and returns the result as a vector
*/
/**************************************************************************/
Vector2D  operator * (const Matrix4x4 &pMtx, const Vector2D &rhs);

/**************************************************************************/
/*!
This function sets the matrix pResult to the identity matrix
*/
/**************************************************************************/
void Mtx44Identity(Matrix4x4 &pResult);

/**************************************************************************/
/*!
This function creates a translation matrix from x & y
and saves it in pResult
*/
/**************************************************************************/
void Mtx44Translate(Matrix4x4 &pResult, float x, float y);

/**************************************************************************/
/*!
This function creates a scaling matrix from x & y
and saves it in pResult
*/
/**************************************************************************/
void Mtx44Scale(Matrix4x4 &pResult, float x, float y);

/**************************************************************************/
/*!
This matrix creates a rotation matrix from "angle" whose value
is in radian. Save the resultant matrix in pResult.
*/
/**************************************************************************/
void Mtx44RotRad(Matrix4x4 &pResult, float angle);

/**************************************************************************/
/*!
This matrix creates a rotation matrix from "angle" whose value
is in degree. Save the resultant matrix in pResult.
*/
/**************************************************************************/
void Mtx44RotDeg(Matrix4x4 &pResult, float angle);

/**************************************************************************/
/*!
This functions calculated the transpose matrix of pMtx
and saves it in pResult
*/
/**************************************************************************/
void Mtx44Transpose(Matrix4x4 &pResult, const Matrix4x4 &pMtx);

/**************************************************************************/
/*!
This function calculates the inverse matrix of pMtx and saves the
result in pResult. If the matrix inversion fails, pResult
would be set to NULL.
*/
/**************************************************************************/
void Mtx44Inverse(Matrix4x4 *pResult, float *determinant, const Matrix4x4 &pMtx);

#endif // CS230_Matrix4X4_H_