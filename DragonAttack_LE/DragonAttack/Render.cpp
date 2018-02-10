/* Start Header ************************************************************************/
/*!
\file       Render.cpp
\author     Jacob Lim
\par email: jacob.lim\@digipen.edu
\brief
  Render class member functions defined here.

Copyright (C) 20xx DigiPen Institute of Technology.
Reproduction or disclosure of this file or its contents
without the prior written consent of DigiPen Institute of
Technology is prohibited.
*/
/* End Header **************************************************************************/

#include "AEEngine.h"
#include "Sprite.h"
#include "Transform.h"
#include "Render.h"

//Renders the given object
void Render::Render_Obj(const Sprite& Obj, const Transform &matrix) const
{
  AEGfxSetRenderMode(Obj.Get_rm());
  AEGfxSetTransform(matrix.GetMatrix().m);
  AEGfxTextureSet(Obj.Get_tex(), Obj.Get_posU(), Obj.Get_posV());
  AEGfxSetTintColor(R, G, B, A);
  AEGfxSetTransparency(Transparency);
  AEGfxSetBlendMode(BM);
  AEGfxMeshDraw(Obj.Get_mesh(), AE_GFX_MDM_TRIANGLES);
}
//Default Constructor sets all to 1
Render::Render()
  :BM{ AE_GFX_BM_NONE }, R{ 1.0f }, G{ 1.0f }, B{ 1.0f }, 
  A{ 1.0f }, Transparency{ 1.0f }
  {}
//This Constructor sets all to user specified values
Render::Render(const float &Red, const float &Green, const float &Blue,
               const float &Alpha, const float &Trans, AEGfxBlendMode BlendMode)
  :BM{ BlendMode }, R{ Red }, G{ Green }, B{ Blue },
  A{ Alpha }, Transparency{ Trans }
  {}
