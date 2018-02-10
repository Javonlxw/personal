/* Start Header ************************************************************************/
/*!
\file       Camera.h
\author     Jacob Lim
\par email: jacob.lim\@digipen.edu
\brief
  Camera functions declaration.

Copyright (C) 20xx DigiPen Institute of Technology.
Reproduction or disclosure of this file or its contents
without the prior written consent of DigiPen Institute of
Technology is prohibited.
*/
/* End Header **************************************************************************/

#pragma once
#include "Transform.h"

/*************************************************************************
Description:
  Make the camera follow an object. 

Parameters:
  PositionX : Set the camera's X position.
  PositionY : Set the camera's Y position.
  OffsetX   : Camera offset in the X direction.
  OffsetY   : Camera offset in the Y direction.

Return:
  void
*************************************************************************/
void CamFollow(float PositionX = 0.0f, float PositionY = 0.0f, float OffsetX = 0.0f, float OffsetY = 0.0f);

/*************************************************************************
Description:
  Make the camera follow an object.

Parameters:
  PosMatrix : Transform matrix to follow
  OffsetX   : Camera offset in the X direction.
  OffsetY   : Camera offset in the Y direction.

Return:
  void
*************************************************************************/
void CamFollow(Transform const &PosMatrix, float OffsetX, float OffsetY, float dir);

/*************************************************************************
Description:
  Set camera position with no offset. 

Parameters:
  PositionX : Set the camera's X position.
  PositionY : Set the camera's Y position.

Return:
  void
*************************************************************************/
void CamStatic(float PositionX = 0.0f, float PositionY = 0.0f);
