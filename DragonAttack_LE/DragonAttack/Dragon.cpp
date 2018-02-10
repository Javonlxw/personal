/* Start Header ************************************************************************/
/*!
\file       Dragon.cpp
\author     Jacob Lim
\par email: jacob.lim\@digipen.edu
\brief
  Dragon class member function definitions here.

Copyright (C) 20xx DigiPen Institute of Technology.
Reproduction or disclosure of this file or its contents
without the prior written consent of DigiPen Institute of
Technology is prohibited.
*/
/* End Header **************************************************************************/

#include "Dragon.h"
#include "Camera.h"
#include "Input_Handler.h"

void Dragon::ApplyPowerUP()
{

}

void Dragon::Input()
{
  const std::vector <int> &Input = Input::Get_User_Input();

  for (unsigned i = 0; i < Input.size(); ++i)
  {
    switch (Input[i])
    {
      case Input::go_left:
        Dir.L = true;
        break;
      case Input::go_right:
        Dir.R = true;
        break;
      case Input::jump_up:
        if (PosY == Start_Pos_Y)
          Dir.UP = true;
        break;
      case Input::fire:
        Attack = true;
        break;
    }
  }
  Input::ClearBuffer();
  //if (AEInputCheckCurr(AEVK_A))
  //{
  //  Dir.L = true;
  //}
  //else if (AEInputCheckCurr(AEVK_D))
  //{
  //  Dir.R = true;
  //}
  //if (AEInputCheckCurr(AEVK_SPACE))
  //  if(PosY == Start_Pos_Y)
  //    Dir.UP = true;
  //if (AEInputCheckTriggered(AEVK_RETURN))
  //  Attack = true;
}

void Dragon::Update()
{
  //Call for player input
  Input();
  //Update position based on direction
  if (Dir.L)
  {
    PosX -= Speed;
    Facing = -1.0f;
    Transform_.SetScale(Facing, 1.0f);
  }
  else if (Dir.R)
  {
    PosX += Speed;
    Facing = 1.0f;
    Transform_.SetScale(Facing, 1.0f);
  }
  if (Dir.UP)
  {
    PosY += Speed * Jump_Mult;
    Air_Dist += Speed * Jump_Mult;
  }
  PosY -= Gravity;
  if (PosY <= Start_Pos_Y)
    PosY = Start_Pos_Y;
  //Update position of player
  Transform_.SetTranslate(PosX, PosY);
  Transform_.Concat();
  Collision_.Update_Col_Pos(PosX - 100.0f, PosY - 100.0f, PosX + 100.0f, PosY + 100.0f);
  //Check if attack has been made
  if (Attack)
  {
    for (int i = 0; i < Bullet_Buffer; ++i)
      if (Fireball[i].IsActive())
      {
        if (Fireball[i].GetDist() <= Bullet_Interval)
          break;
      }
      else
      {
        Fireball[i].SetActive(true);
        if (Facing < 0.0f)
          Fireball[i].SetDir(false);
        else
          Fireball[i].SetDir(true);
        break;
      }
  }
  //Check for any active fireballs
  for (int i = 0; i < Bullet_Buffer; ++i)
  {
    Fireball[i].Pos(PosX, PosY);
    Fireball[i].Pos();
    Fireball[i].Update();
    Fireball[i].Collision_.Update_Col_Pos(Fireball[i].PosX, Fireball[i].PosY);
  }
  //Check for distance limit
  for (int i = 0; i < Bullet_Buffer; ++i)
    if (Fireball[i].IsActive())
    {
      if (Fireball[i].GetDist() >= Bullet_Death)
      {
        Fireball[i].ResetDist();
        Fireball[i].SetActive(false);
      }
    }
    else
      continue;
  //Reset direction and attack
  Attack = false;
  Dir.L = Dir.R = false;
  if (Air_Dist >= Jump_Height)
  {
    Dir.UP = false;
    Air_Dist = 0.0f;
  }
}

void Dragon::Render()
{
  for (int i = 0; i < Bullet_Buffer; ++i)
    Fireball[i].Render();
  Mfireball.Render();
  GameObject::Render();
}
