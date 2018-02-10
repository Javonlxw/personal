/* Start Header ************************************************************************/
/*!
\file       Sprite.cpp
\author     Jacob Lim
\par email: jacob.lim\@digipen.edu
\brief
  Sprite class member functions are defined here.

Copyright (C) 20xx DigiPen Institute of Technology.
Reproduction or disclosure of this file or its contents
without the prior written consent of DigiPen Institute of
Technology is prohibited.
*/
/* End Header **************************************************************************/

#include "Sprite.h"
#include "Transform.h"
#include "AEEngine.h"

//Default constructor; Sets everything to zero
Sprite::Sprite()
  :Mesh{ nullptr }, RM{ AE_GFX_RM_COLOR }, Tex{ nullptr },
  U{ 0.0f }, V{ 0.0f }, Width{ 0.0f }, Height{ 0.0f }
  {}
//Move constructor
Sprite::Sprite(Sprite&& Move_Object)
  : Mesh{ Move_Object.Mesh }, RM{ Move_Object.RM }, Tex{ Move_Object.Tex }, 
  U{ Move_Object.U }, V{ Move_Object.V }, 
  Width{ Move_Object.Width }, Height{ Move_Object.Height }
{
  //Repair Move_Object
  Move_Object.Mesh = nullptr;
  Move_Object.Tex = nullptr;
}
//Move assignment operator
Sprite& Sprite::Sprite::operator=(Sprite&& Move_Object)
{
  if (this != &Move_Object)
  {
    Mesh = Move_Object.Mesh;
    Tex = Move_Object.Tex;
    RM = Move_Object.RM;
    U = Move_Object.U;
    V = Move_Object.V;
    Width = Move_Object.Width;
    Height = Move_Object.Height;
    Move_Object.Mesh = nullptr;
    Move_Object.Tex = nullptr;
  }
  return *this;
}
//Construct with a mesh, but no texture
Sprite::Sprite(AEGfxVertexList * mesh, const float &ObjectW, const float &ObjectH)
  :Mesh{ mesh }, RM{ AE_GFX_RM_COLOR }, Tex{ nullptr },
  U{ 0.0f }, V{ 0.0f }, Width{ ObjectW }, Height{ ObjectH }
  {}
//Construct with a mesh and texture
Sprite::Sprite(AEGfxVertexList * mesh, const char* TexFile, 
               float ObjectW, const float &ObjectH)
  :Mesh{ mesh }, RM{ AE_GFX_RM_TEXTURE }, Tex{ AEGfxTextureLoad(TexFile) }, Draw{},
  U{ 0.0f }, V{ 0.0f }, Width{ ObjectW }, Height{ ObjectH }
  {}
//Free all resources
Sprite::~Sprite()
{
  if (Tex != nullptr)
    AEGfxTextureUnload(Tex);
  if (Mesh != nullptr)
    AEGfxMeshFree(Mesh);
}
//This function renders an object
void Sprite::Render_Object(const Transform &matrix) const
{
  Draw.Render_Obj(*this, matrix);
}
//Set texture position
Sprite& Sprite::SetTexPos(const float &posU, const float &posV)
{
  U = posU;
  V = posV;
  return *this;
}
//Set render parameters
Sprite& Sprite::SetRGB(const float &Red, const float &Green, const float &Blue)
{
  Draw.SetRGB(Red, Green, Blue);
  return *this;
}
//Set Alpha, Transparency, Blend Mode
Sprite& Sprite::SetAlphaTransBM(const float &Alpha, const float &Trans,
                                AEGfxBlendMode BlendMode)
{
  Draw.SetAlpha(Alpha);
  Draw.SetTransparency(Trans);
  Draw.SetBlendMode(BlendMode);
  return *this;
}
