/* Start Header ************************************************************************/
/*!
\file       GameObject.cpp
\author     Jacob Lim
\par email: jacob.lim\@digipen.edu
\brief
  Game Object member functions defined here.

Copyright (C) 20xx DigiPen Institute of Technology.
Reproduction or disclosure of this file or its contents
without the prior written consent of DigiPen Institute of
Technology is prohibited.
*/
/* End Header **************************************************************************/

#include "AEEngine.h"   //AEGfxVertexList
#include "GameObject.h" //Class declaration
#include <utility>      //move

//Renders the object
void GameObject::Render()
{
  if(active)
    Sprite_.Render_Object(Transform_);
}

//Sets the object to active or inactive
void GameObject::SetActive(bool status)
{
  active = status;
}

//Default constructor; Initializes to zero
GameObject::GameObject()
  :active{ false }, Transform_{}, Sprite_{}
{}

//Constructs sprite with a mesh, but no texture
GameObject::GameObject(AEGfxVertexList *mesh, const float &ObjectW, const float &ObjectH)
  : active{ false }, Transform_{}, Sprite_{ mesh, ObjectW, ObjectH }
{}

//Constructs sprite with a mesh and texture
GameObject::GameObject(AEGfxVertexList *mesh, const char* TexFile,
                       float ObjectW, const float &ObjectH)
  : active{ false }, Transform_{}, Sprite_{ mesh, TexFile, ObjectW, ObjectH }
{}

//Move constructor for sprite
GameObject::GameObject(Sprite&& t_Sprite, Col_Comp &&t_Col)
  : active{ false }, Transform_{}, Sprite_{ std::move(t_Sprite) }, Collision_{ std::move(t_Col) }
{}
