/* Start Header ************************************************************************/
/*!
\file       Render.h
\author     Jacob Lim
\par email: jacob.lim\@digipen.edu
\brief
  Render class declared here.

Copyright (C) 20xx DigiPen Institute of Technology.
Reproduction or disclosure of this file or its contents
without the prior written consent of DigiPen Institute of
Technology is prohibited.
*/
/* End Header **************************************************************************/

#pragma once
#include "AEEngine.h"

class Sprite;    //Object foward declaration
class Transform; //Transform foward declaration

/****************************************************************************************
Description: 
  Renders an object.
Constructors:
  Default constructor sets all values to 1.
  Non default constructor sets values by user.
Deconstructor:
  None
****************************************************************************************/
class Render {

public:
  //Set tint RGB
  Render &SetRGB(const float &Red, const float &Green, const float &Blue)
  { R = Red; G = Green; B = Blue; return *this; }
  //Set tint Alpha
  Render &SetAlpha(const float &Alpha) { A = Alpha; return *this; }
  //Set transparency value for rendering
  Render &SetTransparency(const float &Trans) { Transparency = Trans; return *this; }
  //Set blend mode for rendering
  Render &SetBlendMode(AEGfxBlendMode BlendMode) { BM = BlendMode; return *this; }
  //Renders the object
  void Render_Obj(const Sprite &Obj, const Transform &matrix) const;
  //Default Constructor sets all to 1
  Render();
  //This Constructor sets all to user specified values
  Render(const float &Red, const float &Green, const float &Blue, 
         const float &Alpha, const float &Trans, AEGfxBlendMode BlendMode);

private:
  AEGfxBlendMode BM;     //Blend Mode
  float R;               //Red tint value
  float G;               //Green tint value
  float B;               //Blue tint value
  float A;               //Alpha tint value
  float Transparency;    //Transparency value
};
