/* Start Header ************************************************************************/
/*!
\file       GameObject.h
\author     Jacob Lim
\par email: jacob.lim\@digipen.edu
\brief
  Game Object class declared here.

Copyright (C) 20xx DigiPen Institute of Technology.
Reproduction or disclosure of this file or its contents
without the prior written consent of DigiPen Institute of
Technology is prohibited.
*/
/* End Header **************************************************************************/

#pragma once
#include "AEEngine.h"
#include "Transform.h"
#include "Collision.h"
#include "Sprite.h"

/****************************************************************************************
Description:
  Creates a Game Object
Constructors:
  Game Object can be constructed with a mesh or with mesh and texture or none
Deconstructor:
  None
****************************************************************************************/
class GameObject {

public:
  virtual void Update() = 0; //Overwrite with derived class Update
  //Renders the object
  void Render();
  //Sets the object to active or inactive
  void SetActive(bool status);
  //Returns Active or Inactive status
  bool IsActive() const { return active; }
  //Sets object position
  void SetPos(const float &X, const float &Y) { PosX = X; PosY = Y; }
  //Set Velocity
  GameObject &SetVelocity(const AEVec2 &t_vec) { Velocity = t_vec; return *this; }
  //Get Velocity
  AEVec2 GetVelocity() const { return Velocity; }
  //Default constructor; Initializes to zero
  GameObject();
  //Constructs sprite with a mesh, but no texture
  GameObject(AEGfxVertexList *mesh, const float &ObjectW, const float &ObjectH);
  //Constructs sprite with a mesh and texture
  GameObject(AEGfxVertexList * mesh, const char* TexFile,
             float ObjectW, const float &ObjectH);
  //Move constructor for sprite
  GameObject(Sprite&& t_Sprite, Col_Comp &&t_Col);
  Transform Transform_; //Transform holds the object's position, scale and rotation
  Sprite Sprite_;       //Sprite handles mesh, texture and rendering
  Col_Comp Collision_;  //Collision handles object collision
  float PosX{ 0.0f };   //Position X
  float PosY{ 0.0f };   //Position Y

private:
  bool active; //Active or Inactive
  AEVec2 Velocity{ 0.0f, 0.0f };
};
