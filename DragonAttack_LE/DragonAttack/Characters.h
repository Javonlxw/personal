/* Start Header ************************************************************************/
/*!
\file       Characters.h
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

#pragma once
#include <utility>
#include "GameObject.h"
#include "Transform.h"
#include "Create_Object.h"

enum Direction
{
	LEFT = 0,
	RIGHT,
	TOP,
	DOWN,
	STAY
};

// A class that creates a basic character
class Characters : public GameObject
{
	private:
		int hp;      // hp the character have
		boolean moving;  // whether its moving or not
		Direction direction;   // moving in left or right, up or down?

	public:
		Characters(Sprite&& t_sprite, const int HP, Col_Comp &&t_col); // default constructor

		int Get_HP() const;         // get HP of mob
		void Set_HP(int hp);       // set the hp for mob
		void Decrease_HP(const char damage = 1);    // decrease the hp of the character

		boolean Is_Moving() const;     // is character moving
		void Set_Moving(boolean moving); // make a character stop/move

		Direction Get_Direction() const; // get the direction of mob
		void Set_Direction(Direction direction);// set the direction of mob

		bool collision();
};
