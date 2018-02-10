/* Start Header ************************************************************************/
/*!
\file       Sprite.h
\author     Jacob Lim
\par email: jacob.lim\@digipen.edu
\brief
  Sprite class declared here.

Copyright (C) 20xx DigiPen Institute of Technology.
Reproduction or disclosure of this file or its contents
without the prior written consent of DigiPen Institute of
Technology is prohibited.
*/
/* End Header **************************************************************************/

#pragma once
#include "AEEngine.h"
#include "Transform.h"
#include "Render.h"

const float FLOOR_WIDTH { 640.0f }; //The width of the floor
const float FLOOR_HEIGHT{ 100.0f }; //The height of the floor
const float PLAT_WIDTH  { 64.0f  }; //The width of one platform
const float PLAT_HEIGHT { 36.0f  }; //The height of one platform
const float BG_WIDTH    { 640.0f }; //The width of one background
const float BG_HEIGHT   { 360.0f }; //The height of one background

/****************************************************************************************
Description: 
  Creates an object.
Constructors:
  Object can be constructed with a mesh or with mesh and texture or none
Deconstructor:
  Object will free its resources when done
****************************************************************************************/
class Sprite {

public:
  //Set texture position
  Sprite &SetTexPos(const float &posU, const float &posV);
  //Set RGB
  Sprite &SetRGB(const float &Red = 1.0f, const float &Green = 1.0f,
                 const float &Blue = 1.0f);
  //Set Alpha, Transparency, Blend Mode
  Sprite &SetAlphaTransBM(const float &Alpha, const float &Trans,
                          AEGfxBlendMode BlendMode);
  //Renders the object
  void Render_Object(const Transform &matrix) const;
  float Get_posU() const { return U; }
  float Get_posV() const { return V; }
  float Get_Width() const { return Width; }
  float Get_Height() const { return Height; }
  AEGfxVertexList* Get_mesh() const { return Mesh; }
  AEGfxRenderMode Get_rm() const { return RM; }
  AEGfxTexture* Get_tex() const { return Tex; }
  Sprite &Get_this() { return *this; }
  //Default constructor; Sets everything to zero
  Sprite();
  //Move constructor
  Sprite(Sprite&& Move_Object);
  //Move assignment operator
  Sprite& Sprite::operator=(Sprite&& Move_Object);
  //Construct with a mesh, but no texture
  Sprite(AEGfxVertexList * mesh, const float &ObjectW, const float &ObjectH);
  //Construct with a mesh and texture
  Sprite(AEGfxVertexList * mesh, const char* TexFile,
         float ObjectW, const float &ObjectH);
  //Free all resources
  ~Sprite();

private:
  AEGfxVertexList *Mesh; //A pointer to object's mesh
  AEGfxRenderMode RM;    //Render Mode
  AEGfxTexture *Tex;     //Texture
  Render Draw;           //Draw component
  float U;               //Texture U
  float V;               //Texture V
  float Width;           //Width of object
  float Height;          //Height of object
};
