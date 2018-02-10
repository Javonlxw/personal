/* Start Header ************************************************************************/
/*!
\file       Projectile.h
\author     Jacob Lim
\par email: jacob.lim\@digipen.edu
\brief
  Projectile class declared here.

Copyright (C) 20xx DigiPen Institute of Technology.
Reproduction or disclosure of this file or its contents
without the prior written consent of DigiPen Institute of
Technology is prohibited.
*/
/* End Header **************************************************************************/

#pragma once
#include "GameObject.h"
#include "f_sqrt.h"
#include <utility>

class Projectile : public GameObject {

public:
  //Updates the projectile
  void Update() 
  {
    //Update distance travelled and resultant matrix
    if (IsActive())
    {
      Transform_.SetTranslate(PosX, PosY);
      Collision_.Update_Col_Pos(PosX, PosY);
      //Pythagoras distance equation
      AddDist(
       f_sqrt( (GetVelocity().x) * 0.016f * (GetVelocity().x * 0.016f)
                                +
              (GetVelocity().y) * 0.016f * (GetVelocity().y) * 0.016f)
      );
      Transform_.Concat();
    }
  }
  //Set the position of the projectile
  void Pos()
  {
    //Check for active
    if (IsActive())
    {
      if (GetDir())
      {
        PosX += GetVelocity().x * 0.016f;
        PosY += GetVelocity().y * 0.016f;
        Transform_.SetRotation(90.0f);
      }
      else
      {
        PosX -= GetVelocity().x;
        PosY -= GetVelocity().y;
        Transform_.SetRotation(-90.0f);
      }
    }
  }
  //Updates position of projectile if not active
  void Pos(const float &X, const float &Y)
  {
    if (!IsActive())
    {
      PosX = X;
      PosY = Y;
    }
  }
  //Distance travelled by the projectile
  void AddDist(const float &dist) { Distance += dist; }
  //Returns Distance travelled by the projectile
  float GetDist() const { return Distance; }
  //Resets the distance travelled
  void ResetDist() { Distance = 0.0f; }
  //Sets the direction of the projectile
  void SetDir(const bool &dir) { Direction = dir; }
  //Returns the direction of the projectile
  bool GetDir() const { return Direction; }
  //Default contructor does nothing
  Projectile() = default;
  //Initialize projectile with mesh and collision shape
  Projectile(Sprite &&t_sprite, Col_Comp &&t_col)
    :GameObject{ std::move(t_sprite), std::move(t_col) }, 
    Distance{ 0.0f }, Direction{ true }
  {}
private:
  bool Direction; //TRUE = Right; FALSE = Left
  float Distance; //Distance travelled
};
