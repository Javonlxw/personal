/* Start Header ************************************************************************/
/*!
\file       Create_Object.cpp
\author     Jacob Lim
\par email: jacob.lim\@digipen.edu
\brief
  Mesh creation functions are defined here.

Copyright (C) 20xx DigiPen Institute of Technology.
Reproduction or disclosure of this file or its contents
without the prior written consent of DigiPen Institute of
Technology is prohibited.
*/
/* End Header **************************************************************************/

#include "AEEngine.h"
#include "Create_Object.h"
#include "Sprite.h"

//This function creates a custom square
AEGfxVertexList* CreateSquare(float size, float scaleU, float scaleV, unsigned int color)
{
  AEGfxMeshStart();
  AEGfxTriAdd(
    -size, -size, color, 0.0f, 1.0f * scaleV,
    size, -size, color, 1.0f * scaleU, 1.0f * scaleU,
    -size, size, color, 0.0f, 0.0f);
  AEGfxTriAdd(
    size, -size, color, 1.0f * scaleU, 1.0f * scaleV,
    size, size, color, 1.0f * scaleU, 0.0f,
    -size, size, color, 0.0f, 0.0f);

  return AEGfxMeshEnd();
}

//This function creates a custom rectangle
AEGfxVertexList* CreateRectangle(float width, float height, float ScaleU, float ScaleV, unsigned int color)
{
  AEGfxMeshStart();
  AEGfxTriAdd(
    -width, -height, color, 0.0f, 1.0f,
    width, -height, color, (width / height) * ScaleU, 1.0f,
    -width, height, color, 0.0f, -1.0f * ScaleV + 1.0f);
  AEGfxTriAdd(
    width, -height, color, 1.0f, 1.0f * ScaleV,
    width, height, color, 1.0f, 0.0f,
    -width, height, color, -(width / height) * ScaleU + 1.0f, 0.0f);

  return AEGfxMeshEnd();
}

//This function returns a sprite object with a square mesh
Sprite S_CreateSquare(float size, float scaleU, float scaleV, const char* texFile, unsigned int color)
{
  if (texFile)
  {
    Sprite tmp_Square(CreateSquare(size, scaleU, scaleV, color), texFile, size, size);
    return tmp_Square;
  }
  else
  {
    Sprite tmp_Square(CreateSquare(size, scaleU, scaleV, color), size, size);
    return tmp_Square;
  }
}

//This function returns a sprite object with a rectangle mesh
Sprite S_CreateRectangle(float width, float height, float ScaleU, float ScaleV, const char* texFile, unsigned int color)
{
  if (texFile)
  {
    Sprite tmp_Rectangle(CreateRectangle(width, height, ScaleU, ScaleV, color), texFile, width, height);
    return tmp_Rectangle;
  }
  else
  {
    Sprite tmp_Rectangle(CreateRectangle(width, height, ScaleU, ScaleV, color), width, height);
    return tmp_Rectangle;
  }
}

//This function returns a sprite object with a platform mesh
Sprite CreatePlatform(float MultW, float MultH, float ScaleU, float ScaleV, const char* texFile, unsigned int color)
{
  if (texFile)
  {
    Sprite tmp_Rectangle(
      CreateRectangle(PLAT_WIDTH * MultW, PLAT_HEIGHT * MultH, ScaleU, ScaleV, color), 
      texFile, PLAT_WIDTH * MultW, PLAT_HEIGHT * MultH);
    return tmp_Rectangle;
  }
  else
  {
    Sprite tmp_Rectangle(CreateRectangle(
      PLAT_WIDTH * MultW, PLAT_HEIGHT * MultH, ScaleU, ScaleV, color), 
      PLAT_WIDTH * MultW, PLAT_HEIGHT * MultH);
    return tmp_Rectangle;
  }
}

//This function returns a sprite object with a floor mesh
Sprite CreateFloor(float MultW, float ScaleU, float ScaleV, const char* texFile, unsigned int color)
{
  if (texFile)
  {
    Sprite tmp_Rectangle(
      CreateRectangle(FLOOR_WIDTH * MultW, FLOOR_HEIGHT, ScaleU, ScaleV, color),
      texFile, FLOOR_WIDTH * MultW, FLOOR_HEIGHT);
    return tmp_Rectangle;
  }
  else
  {
    Sprite tmp_Rectangle(CreateRectangle(
      FLOOR_WIDTH * MultW, FLOOR_HEIGHT, ScaleU, ScaleV, color),
      FLOOR_WIDTH * MultW, FLOOR_HEIGHT);
    return tmp_Rectangle;
  }
}

//This function returns a sprite object with a BG mesh
Sprite CreateBG(float MultW, float ScaleU, float ScaleV, const char* texFile, unsigned int color)
{
  if (texFile)
  {
    Sprite tmp_Rectangle(
      CreateRectangle(BG_WIDTH * MultW, BG_HEIGHT, ScaleU, ScaleV, color),
      texFile, BG_WIDTH * MultW, BG_HEIGHT);
    return tmp_Rectangle;
  }
  else
  {
    Sprite tmp_Rectangle(CreateRectangle(
      BG_WIDTH * MultW, BG_HEIGHT, ScaleU, ScaleV, color),
      BG_WIDTH * MultW, BG_HEIGHT);
    return tmp_Rectangle;
  }
}
