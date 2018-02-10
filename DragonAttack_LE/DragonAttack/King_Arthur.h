/* Start Header ************************************************************************/
/*!
\file       King_Arthur.h
\author     William Yoong
\par email: william.yoong\@digipen.edu
\brief
This file creates, set the behaviour and rules for the last boss king arthur

Copyright (C) 2018 DigiPen Institute of Technology.
Reproduction or disclosure of this file or its contents
without the prior written consent of DigiPen Institute of
Technology is prohibited.
*/
/* End Header **************************************************************************/
#pragma once

#include "Characters.h"
#include "Dragon.h"

class King_Arthur : public Characters
{
private:
	enum KA_MoveSet // for switch case to determine which attack gets called
	{
		SINGLE_SLASH = 0,
		TRIPLE_SLASH,
		UNIQUE_MECHANIC
	};
	
	void Idle(float dt);

	void Moving(const Dragon &d, float dt);

	void AvoidingObstacle();

	void Attack(const Dragon &d, float dt);

	void Init(void);

	void King_Arthur_Phase2(void); // changes mechanics for phase 2

	void Jump_Attack(void); // jump attack

	void Single_Slash(void); // single slash + wave

	void Triple_Slash(void); // single slash 3 times

	void Heal_and_Spawn(void); // phase 2, heal and spawm mobs

	void(King_Arthur::*ka_attacks[UNIQUE_MECHANIC + 1])(void);

	bool phase2;

public:
	/******************************************************************************/
	/*!
	\fn         King Arthuer

	\brief      Constructor for King Arthur

	\param      Sprite&& t_sprite
	- an rvalue of the sprite to be moved

	\retval     void
	No return.
	*/
	/******************************************************************************/
	King_Arthur();

	~King_Arthur(); // destructor
	void Update(float dt, const Dragon &player);
	void Update() override {}; // ignored

};

