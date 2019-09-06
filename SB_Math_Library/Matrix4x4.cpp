#include "Matrix4x4.h"
#include <cmath>
#include <glm/gtc/type_ptr.hpp>
#define EPSILON		0.0001f
#define PI			3.14159265358f

/*****************************************************************************/
/*!
Constructor
*/
/*****************************************************************************/
Matrix4x4::Matrix4x4(const float *pArr)
{
  //Initiliaze the m[16] array
  for (int i = 0; i < 4; ++i)
    for (int k = 0; k < 4; ++k)
      m[i][k] = *(pArr + i + k);
}
/*****************************************************************************/
/*!
Constructor
*/
/*****************************************************************************/
Matrix4x4::Matrix4x4(float _00, float _01, float _02, float _03,
  float _10, float _11, float _12, float _13,
  float _20, float _21, float _22, float _23,
  float _30, float _31, float _32, float _33)
{
  m[0][0] = _00;
  m[0][1] = _01;
  m[0][2] = _02;
  m[0][3] = _03;
  m[1][0] = _10;
  m[1][1] = _11;
  m[1][2] = _12;
  m[1][3] = _13;
  m[2][0] = _20;
  m[2][1] = _21;
  m[2][2] = _22;
  m[2][3] = _23;
  m[3][0] = _30;
  m[3][1] = _31;
  m[3][2] = _32;
  m[3][3] = _33;
}

Matrix4x4::Matrix4x4(const glm::mat4& rhs)
{
  m[0][0] = rhs[0][0];
  m[0][1] = rhs[0][1];
  m[0][2] = rhs[0][2];
  m[0][3] = rhs[0][3];
  m[1][0] = rhs[1][0];
  m[1][1] = rhs[1][1];
  m[1][2] = rhs[1][2];
  m[1][3] = rhs[1][3];
  m[2][0] = rhs[2][0];
  m[2][1] = rhs[2][1];
  m[2][2] = rhs[2][2];
  m[2][3] = rhs[2][3];
  m[3][0] = rhs[3][0];
  m[3][1] = rhs[3][1];
  m[3][2] = rhs[3][2];
  m[3][3] = rhs[3][3];
}

Matrix4x4::Matrix4x4(const Matrix4x4& rhs)
{
	for (int i = 0; i < 4; ++i)
		for (int k = 0; k < 4; ++k)
			m[i][k] = rhs.m[k][i];
}

Matrix4x4::Matrix4x4(float x)
{
	m[0][0] = x;
	m[0][1] = 0;
	m[0][2] = 0;
	m[0][3] = 0;
	m[1][0] = 0;
	m[1][1] = x;
	m[1][2] = 0;
	m[1][3] = 0;
	m[2][0] = 0;
	m[2][1] = 0;
	m[2][2] = x;
	m[2][3] = 0;
	m[3][0] = 0;
	m[3][1] = 0;
	m[3][2] = 0;
	m[3][3] = x;
}

/*****************************************************************************/
/*!
operator=
*/
/*****************************************************************************/
Matrix4x4& Matrix4x4::operator=(const Matrix4x4 &rhs)
{
  for (int i = 0; i < 4; ++i)
    for (int k = 0; k < 4; ++k)
      m[i][k] = rhs.m[i][k];

  return *this;
}
/*****************************************************************************/
/*!
operator*=
*/
/*****************************************************************************/
Matrix4x4& Matrix4x4::operator*=(const Matrix4x4 &rhs)
{
  // Multiply the first row
  m[0][0] = m[0][0] * rhs.m[0][0] + m[0][1] * rhs.m[1][0] + m[0][2] * rhs.m[2][0] + m[0][3] * rhs.m[3][0];
  m[0][1] = m[0][0] * rhs.m[0][1] + m[0][1] * rhs.m[1][1] + m[0][2] * rhs.m[2][1] + m[0][3] * rhs.m[3][1];
  m[0][2] = m[0][0] * rhs.m[0][2] + m[0][1] * rhs.m[1][2] + m[0][2] * rhs.m[2][2] + m[0][3] * rhs.m[3][2];
  m[0][3] = m[0][0] * rhs.m[0][3] + m[0][1] * rhs.m[1][3] + m[0][2] * rhs.m[2][3] + m[0][3] * rhs.m[3][3];
  // Multiply the second row
  m[1][0] = m[1][0] * rhs.m[0][0] + m[1][1] * rhs.m[1][0] + m[1][2] * rhs.m[2][0] + m[1][3] * rhs.m[3][0];
  m[1][1] = m[1][0] * rhs.m[0][1] + m[1][1] * rhs.m[1][1] + m[1][2] * rhs.m[2][1] + m[1][3] * rhs.m[3][1];
  m[1][2] = m[1][0] * rhs.m[0][2] + m[1][1] * rhs.m[1][2] + m[1][2] * rhs.m[2][2] + m[1][3] * rhs.m[3][2];
  m[1][3] = m[1][0] * rhs.m[0][3] + m[1][1] * rhs.m[1][3] + m[1][2] * rhs.m[2][3] + m[1][3] * rhs.m[3][3];
  // Multiply the third row
  m[2][0] = m[2][0] * rhs.m[0][0] + m[2][1] * rhs.m[1][0] + m[2][2] * rhs.m[2][0] + m[2][3] * rhs.m[3][0];
  m[2][1] = m[2][0] * rhs.m[0][1] + m[2][1] * rhs.m[1][1] + m[2][2] * rhs.m[2][1] + m[2][3] * rhs.m[3][1];
  m[2][2] = m[2][0] * rhs.m[0][2] + m[2][1] * rhs.m[1][2] + m[2][2] * rhs.m[2][2] + m[2][3] * rhs.m[3][2];
  m[2][3] = m[2][0] * rhs.m[0][3] + m[2][1] * rhs.m[1][3] + m[2][2] * rhs.m[2][3] + m[2][3] * rhs.m[3][3];
  // Multiply the fourth row
  m[3][0] = m[3][0] * rhs.m[0][0] + m[3][1] * rhs.m[1][0] + m[3][2] * rhs.m[2][0] + m[3][3] * rhs.m[3][0];
  m[3][1] = m[3][0] * rhs.m[0][1] + m[3][1] * rhs.m[1][1] + m[3][2] * rhs.m[2][1] + m[3][3] * rhs.m[3][1];
  m[3][2] = m[3][0] * rhs.m[0][2] + m[3][1] * rhs.m[1][2] + m[3][2] * rhs.m[2][2] + m[3][3] * rhs.m[3][2];
  m[3][3] = m[3][0] * rhs.m[0][3] + m[3][1] * rhs.m[1][3] + m[3][2] * rhs.m[2][3] + m[3][3] * rhs.m[3][3];

  return *this;
}
/*****************************************************************************/
/*!
operator*
*/
/*****************************************************************************/
Matrix4x4 operator*(const Matrix4x4 &lhs, const Matrix4x4 &rhs)
{
  Matrix4x4 temp; //Create a temp Matrix4x4 to return
                  // Multiply the first row
  temp.m[0][0] = lhs.m[0][0] * rhs.m[0][0] + lhs.m[0][1] * rhs.m[1][0] + lhs.m[0][2] * rhs.m[2][0] + lhs.m[0][3] * rhs.m[3][0];
  temp.m[0][1] = lhs.m[0][0] * rhs.m[0][1] + lhs.m[0][1] * rhs.m[1][1] + lhs.m[0][2] * rhs.m[2][1] + lhs.m[0][3] * rhs.m[3][1];
  temp.m[0][2] = lhs.m[0][0] * rhs.m[0][2] + lhs.m[0][1] * rhs.m[1][2] + lhs.m[0][2] * rhs.m[2][2] + lhs.m[0][3] * rhs.m[3][2];
  temp.m[0][3] = lhs.m[0][0] * rhs.m[0][3] + lhs.m[0][1] * rhs.m[1][3] + lhs.m[0][2] * rhs.m[2][3] + lhs.m[0][3] * rhs.m[3][3];
  // Multiply the second row
  temp.m[1][0] = lhs.m[1][0] * rhs.m[0][0] + lhs.m[1][1] * rhs.m[1][0] + lhs.m[1][2] * rhs.m[2][0] + lhs.m[1][3] * rhs.m[3][0];
  temp.m[1][1] = lhs.m[1][0] * rhs.m[0][1] + lhs.m[1][1] * rhs.m[1][1] + lhs.m[1][2] * rhs.m[2][1] + lhs.m[1][3] * rhs.m[3][1];
  temp.m[1][2] = lhs.m[1][0] * rhs.m[0][2] + lhs.m[1][1] * rhs.m[1][2] + lhs.m[1][2] * rhs.m[2][2] + lhs.m[1][3] * rhs.m[3][2];
  temp.m[1][3] = lhs.m[1][0] * rhs.m[0][3] + lhs.m[1][1] * rhs.m[1][3] + lhs.m[1][2] * rhs.m[2][3] + lhs.m[1][3] * rhs.m[3][3];
  // Multiply the third row
  temp.m[2][0] = lhs.m[2][0] * rhs.m[0][0] + lhs.m[2][1] * rhs.m[1][0] + lhs.m[2][2] * rhs.m[2][0] + lhs.m[2][3] * rhs.m[3][0];
  temp.m[2][1] = lhs.m[2][0] * rhs.m[0][1] + lhs.m[2][1] * rhs.m[1][1] + lhs.m[2][2] * rhs.m[2][1] + lhs.m[2][3] * rhs.m[3][1];
  temp.m[2][2] = lhs.m[2][0] * rhs.m[0][2] + lhs.m[2][1] * rhs.m[1][2] + lhs.m[2][2] * rhs.m[2][2] + lhs.m[2][3] * rhs.m[3][2];
  temp.m[2][3] = lhs.m[2][0] * rhs.m[0][3] + lhs.m[2][1] * rhs.m[1][3] + lhs.m[2][2] * rhs.m[2][3] + lhs.m[2][3] * rhs.m[3][3];
  // Multiply the fourth row
  temp.m[3][0] = lhs.m[3][0] * rhs.m[0][0] + lhs.m[3][1] * rhs.m[1][0] + lhs.m[3][2] * rhs.m[2][0] + lhs.m[3][3] * rhs.m[3][0];
  temp.m[3][1] = lhs.m[3][0] * rhs.m[0][1] + lhs.m[3][1] * rhs.m[1][1] + lhs.m[3][2] * rhs.m[2][1] + lhs.m[3][3] * rhs.m[3][1];
  temp.m[3][2] = lhs.m[3][0] * rhs.m[0][2] + lhs.m[3][1] * rhs.m[1][2] + lhs.m[3][2] * rhs.m[2][2] + lhs.m[3][3] * rhs.m[3][2];
  temp.m[3][3] = lhs.m[3][0] * rhs.m[0][3] + lhs.m[3][1] * rhs.m[1][3] + lhs.m[3][2] * rhs.m[2][3] + lhs.m[3][3] * rhs.m[3][3];

  return temp;
}
/*****************************************************************************/
/*!
operator* //TODO: check vector with 4 by 4
*/
/*****************************************************************************/
Vector2D operator*(const Matrix4x4 &pMtx, const Vector2D &rhs)
{
  Vector2D temp; //Create a temp Vector2D to return
                 //Homogenous 2D
  temp.x = rhs.x * pMtx.m[0][0] + rhs.y * pMtx.m[0][1] + pMtx.m[0][2];
  temp.y = rhs.x * pMtx.m[1][0] + rhs.y * pMtx.m[1][1] + pMtx.m[1][2];

  return temp;
}
/*****************************************************************************/
/*!
Mtx44Identity
*/
/*****************************************************************************/
void Mtx44Identity(Matrix4x4 &pResult)
{
  //Creating Identity matrix
  /*
  1 0 0 0
  0 1 0 0
  0 0 1 0
  0 0 0 1
  */
  pResult.m[0][0] = 1.f;
  pResult.m[0][1] = 0.f;
  pResult.m[0][2] = 0.f;
  pResult.m[0][3] = 0.f;

  pResult.m[1][0] = 0.f;
  pResult.m[1][1] = 1.f;
  pResult.m[1][2] = 0.f;
  pResult.m[1][3] = 0.f;

  pResult.m[2][0] = 0.f;
  pResult.m[2][1] = 0.f;
  pResult.m[2][2] = 1.f;
  pResult.m[2][3] = 0.f;

  pResult.m[3][0] = 0.f;
  pResult.m[3][1] = 0.f;
  pResult.m[3][2] = 0.f;
  pResult.m[3][3] = 1.f;
}
/*****************************************************************************/
/*!
Mtx44Translate
*/
/*****************************************************************************/
void Mtx44Translate(Matrix4x4 &pResult, float x, float y)
{
  pResult.m[0][0] = 1;
  pResult.m[0][1] = 0;
  pResult.m[0][2] = 0; //Translating by x amount
  pResult.m[0][3] = 0;

  pResult.m[1][0] = 0;
  pResult.m[1][1] = 1;
  pResult.m[1][2] = 0; //Translating by y amount
  pResult.m[1][3] = 0;

  pResult.m[2][0] = 0;
  pResult.m[2][1] = 0;
  pResult.m[2][2] = 1;
  pResult.m[2][3] = 0;

  pResult.m[3][0] = x;
  pResult.m[3][1] = y;
  pResult.m[3][2] = 0;
  pResult.m[3][3] = 1;
}
/*****************************************************************************/
/*!
Mtx44Scale
*/
/*****************************************************************************/
void Mtx44Scale(Matrix4x4 &pResult, float x, float y)
{
  pResult.m[0][0] = x; //Scaling by x amount
  pResult.m[0][1] = 0;
  pResult.m[0][2] = 0;
  pResult.m[0][3] = 0;

  pResult.m[1][0] = 0;
  pResult.m[1][1] = y; //Scaling by y amount
  pResult.m[1][2] = 0;
  pResult.m[1][3] = 0;

  pResult.m[2][0] = 0;
  pResult.m[2][1] = 0;
  pResult.m[2][2] = 1; //Scaling by 1
  pResult.m[2][3] = 0;

  pResult.m[3][0] = 0;
  pResult.m[3][1] = 0;
  pResult.m[3][2] = 0;
  pResult.m[3][3] = 1; //Scaling by 1
}
/*****************************************************************************/
/*!
Mtx44RotRad
*/
/*****************************************************************************/
void Mtx44RotRad(Matrix4x4 &pResult, float angle)
{
  // Creating the Rotation matrix
  /*
  cos -sin 0 0
  sin  cos 0 0
  0    0   1 0
  0    0   0 1
  */
  pResult.m[0][0] = cos(angle);
  pResult.m[0][1] = -sin(angle);
  pResult.m[0][2] = 0;
  pResult.m[0][3] = 0;

  pResult.m[1][0] = sin(angle);
  pResult.m[1][1] = cos(angle);
  pResult.m[1][2] = 0;
  pResult.m[1][3] = 0;

  pResult.m[2][0] = 0;
  pResult.m[2][1] = 0;
  pResult.m[2][2] = 1;
  pResult.m[2][3] = 0;

  pResult.m[3][0] = 0;
  pResult.m[3][1] = 0;
  pResult.m[3][2] = 0;
  pResult.m[3][3] = 1;
}
/*****************************************************************************/
/*!
Mtx44RotDeg
*/
/*****************************************************************************/
void Mtx44RotDeg(Matrix4x4 &pResult, float angle)
{
  // Creating the Rotation matrix
  /*
  cos -sin 0 0
  sin  cos 0 0
  0    0   1 0
  0    0   0 1
  */
  float angleRad = angle * PI / 180.f; // Converting it to Radians
  pResult.m[0][0] = cos(angleRad);
  pResult.m[0][1] = sin(angleRad);
  pResult.m[0][2] = 0;
  pResult.m[0][3] = 0;

  pResult.m[1][0] = -sin(angleRad);
  pResult.m[1][1] = cos(angleRad);
  pResult.m[1][2] = 0;
  pResult.m[1][3] = 0;

  pResult.m[2][0] = 0;
  pResult.m[2][1] = 0;
  pResult.m[2][2] = 1;
  pResult.m[2][3] = 0;

  pResult.m[3][0] = 0;
  pResult.m[3][1] = 0;
  pResult.m[3][2] = 0;
  pResult.m[3][3] = 1;
}
/*****************************************************************************/
/*!
Mtx44Transpose
*/
/*****************************************************************************/
void Mtx44Transpose(Matrix4x4 &pResult, const Matrix4x4 &pMtx)
{
  //Swapping the rows and columns
  /*
  a b c t      a d g x
  d e f j -->  b e h y
  g h i k 	 c f i z
  x y z q      t j k q
  */
  pResult.m[0][1] = pMtx.m[1][0];
  pResult.m[1][0] = pMtx.m[0][1];
  pResult.m[0][3] = pMtx.m[3][0];
  pResult.m[3][0] = pMtx.m[0][3];
  pResult.m[0][0] = pMtx.m[0][0];

  pResult.m[2][0] = pMtx.m[0][2];
  pResult.m[0][2] = pMtx.m[2][0];
  pResult.m[2][3] = pMtx.m[3][2];
  pResult.m[3][2] = pMtx.m[2][3];
  pResult.m[1][1] = pMtx.m[1][1];

  pResult.m[2][1] = pMtx.m[1][2];
  pResult.m[1][2] = pMtx.m[2][1];
  pResult.m[2][2] = pMtx.m[2][2];

  pResult.m[3][1] = pMtx.m[1][3];
  pResult.m[1][3] = pMtx.m[3][1];
  pResult.m[3][3] = pMtx.m[3][3];
}
/*****************************************************************************/
/*!
Mtx44Inverse //TODO: check whether is it usable
*/
/*****************************************************************************/
void Mtx44Inverse(Matrix4x4 *pResult, float *determinant, const Matrix4x4 &pMtx)
{
  // Finding the determinant
  float FirstMinor = pMtx.m[0][0] * (pMtx.m[1][1] * pMtx.m[2][2] - pMtx.m[2][1] * pMtx.m[1][2]);
  float SecondMinor = -pMtx.m[0][1] * (pMtx.m[1][0] * pMtx.m[2][2] - pMtx.m[2][0] * pMtx.m[1][2]);
  float ThirdMinor = pMtx.m[0][2] * (pMtx.m[1][0] * pMtx.m[2][1] - pMtx.m[2][0] * pMtx.m[1][1]);
  *determinant = FirstMinor + SecondMinor + ThirdMinor;

  if (*determinant == 0) //if determinant is 0, Mtx is not invertible
    pResult = nullptr;
  else
  {
    Matrix4x4 temp;
    temp = pMtx; //To prevent overlapping/using of the same matrix
                 //Transpose the original matrix
    Mtx44Transpose(*pResult, temp);
    temp = *pResult;
    //Find the determinant of each of the 2x2 minor matrices
    pResult->m[0][0] = temp.m[1][1] * temp.m[2][2] - temp.m[2][1] * temp.m[1][2];
    pResult->m[0][1] = temp.m[1][0] * temp.m[2][2] - temp.m[2][0] * temp.m[1][2];
    pResult->m[0][2] = temp.m[1][0] * temp.m[2][1] - temp.m[2][0] * temp.m[1][1];
    pResult->m[1][0] = temp.m[0][1] * temp.m[2][2] - temp.m[2][1] * temp.m[0][2];
    pResult->m[1][1] = temp.m[0][0] * temp.m[2][2] - temp.m[2][0] * temp.m[0][2];

    pResult->m[1][2] = temp.m[0][0] * temp.m[2][1] - temp.m[2][0] * temp.m[0][1];

    pResult->m[2][0] = temp.m[0][1] * temp.m[1][2] - temp.m[1][1] * temp.m[0][2];
    pResult->m[2][1] = temp.m[0][0] * temp.m[1][2] - temp.m[1][0] * temp.m[0][2];
    pResult->m[2][2] = temp.m[0][0] * temp.m[1][1] - temp.m[1][0] * temp.m[0][1];

    //Create the matrix of cofactors. Hence the inverse
    pResult->m[0][0] /= *determinant;
    pResult->m[0][1] /= -*determinant;
    pResult->m[0][2] /= *determinant;

    pResult->m[1][0] /= -*determinant;
    pResult->m[1][1] /= *determinant;
    pResult->m[1][2] /= -*determinant;

    pResult->m[2][0] /= *determinant;
    pResult->m[2][1] /= -*determinant;
    pResult->m[2][2] /= *determinant;
  }
}

glm::mat4 Matrix4x4::ConvertMtx44toGlmMat4()
{
  glm::mat4 tmp;

  tmp[0][0] = m[0][0];
  tmp[0][1] = m[0][1];
  tmp[0][2] = m[0][2];
  tmp[0][3] = m[0][3];
  tmp[1][0] = m[1][0];
  tmp[1][1] = m[1][1];
  tmp[1][2] = m[1][2];
  tmp[1][3] = m[1][3];
  tmp[2][0] = m[2][0];
  tmp[2][1] = m[2][1];
  tmp[2][2] = m[2][2];
  tmp[2][3] = m[2][3];
  tmp[3][0] = m[3][0];
  tmp[3][1] = m[3][1];
  tmp[3][2] = m[3][2];
  tmp[3][3] = m[3][3];

  return tmp;
}

