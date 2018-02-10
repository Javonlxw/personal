/* Start Header ************************************************************************/
/*!
\file       Transform.h
\author     Jacob Lim
\par email: jacob.lim\@digipen.edu
\brief
  Transform class declared here.

Copyright (C) 20xx DigiPen Institute of Technology.
Reproduction or disclosure of this file or its contents
without the prior written consent of DigiPen Institute of
Technology is prohibited.
*/
/* End Header **************************************************************************/

#pragma once
#include "AEEngine.h"

/****************************************************************************************
Description:
  Creates a transform class.
Constructors:
  Sets all matrices to the Identity Matrix.
Deconstructor:
  None.
****************************************************************************************/
class Transform {
public:
  //Sets the translation matrix and applies it to the matrix
  Transform &SetTranslate(float posX, float posY);
  //Sets the degrees to rotate and applies it to the matrix
  Transform &SetRotation(const float &Deg);
  //Sets the scale matrix and applies it to the matrix
  Transform &SetScale(float scaleX, float scaleY);
  void Concat() { 
    AEMtx33Concat(&Matrix, &S_Matrix, &R_Matrix); 
    AEMtx33Concat(&Matrix, &T_Matrix, &Matrix);
  }
  float GetDeterminant() { return AEMtx33Determinant(&Matrix); }
  float GetDegree() const { return Degree; }
  AEMtx33 GetMatrix() const { return Matrix; }
  AEMtx33 GetTranslateMatrix() const { return T_Matrix; }
  AEMtx33 GetScaleMatrix() const { return T_Matrix; }
  //Sets all matrices to the Identity Matrix
  Transform()
    :Matrix{ 0.0f }, T_Matrix{ 0.0f }, S_Matrix{ 0.0f }, R_Matrix{ 0.0f }, Degree{ 0.0f }
  {
    AEMtx33Identity(&Matrix);
    AEMtx33Identity(&T_Matrix);
    AEMtx33Identity(&R_Matrix);
    AEMtx33Identity(&S_Matrix);
  }
private:
  AEMtx33 Matrix;   //Resultant matrix
  AEMtx33 T_Matrix; //Translation matrix
  AEMtx33 R_Matrix; //Scale matrix
  AEMtx33 S_Matrix; //Scale matrix
  float Degree;     //Rotation in degrees
};
