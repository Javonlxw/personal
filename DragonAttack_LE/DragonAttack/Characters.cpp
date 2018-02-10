/* Start Header ************************************************************************/
/*!
\file       Characters.cpp
\author     William Yoong
\par email: william.yoong\@digipen.edu
\brief
Character class, parent class of the characters to be created

Copyright (C) 20xx DigiPen Institute of Technology.
Reproduction or disclosure of this file or its contents
without the prior written consent of DigiPen Institute of
Technology is prohibited.
*/
/* End Header **************************************************************************/

#include "Characters.h"
#include <utility>


Characters::Characters(Sprite&& t_sprite, const int HP, Col_Comp &&t_col)
	: GameObject{ std::move(t_sprite), std::move(t_col) }, hp{ HP }, moving{ false }, 
	  direction{ RIGHT }
{
}

boolean Characters::Is_Moving() const
{
	return moving;
}

void Characters::Set_Moving(boolean movement)
{
	moving = movement;
}

void Characters::Set_HP(int hit_points)
{
	hp = hit_points;
}

void Characters::Decrease_HP(const char damage)
{
	hp -= (int)damage;
}

int Characters::Get_HP() const
{
	return hp;
}

Direction Characters::Get_Direction(void) const
{
	return direction;
}

void Characters::Set_Direction(Direction d)
{
	direction = d;
}

bool Characters::collision()
{
	return false;
}