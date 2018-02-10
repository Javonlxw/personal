/* Start Header ************************************************************************/
/*!
\file       Transform.cpp
\author     Jacob Lim
\par email: jacob.lim\@digipen.edu
\brief
  Transform class functions are defined here.

Copyright (C) 20xx DigiPen Institute of Technology.
Reproduction or disclosure of this file or its contents
without the prior written consent of DigiPen Institute of
Technology is prohibited.
*/
/* End Header **************************************************************************/

#include "AEEngine.h"
#include "Transform.h"

//Sets the translation matrix and applies it to the matrix
Transform& Transform::SetTranslate(float posX, float posY)
{
  AEMtx33Trans(&T_Matrix, posX, posY);
  return *this;
}
//Sets the degrees to rotate and applies it to the matrix
Transform& Transform::SetRotation(const float &Deg)
{
  Degree = Deg;
  AEMtx33RotDeg(&R_Matrix, Degree);
  return *this;
}
//Sets the scale matrix and applies it to the matrix
Transform& Transform::SetScale(float scaleX, float scaleY)
{
  AEMtx33Scale(&S_Matrix, scaleX, scaleY);
  return *this;
}