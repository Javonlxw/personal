/* Start Header ************************************************************************/
/*!
\file       Lancelot.h
\author     William Yoong
\par email: william.yoong\@digipen.edu
\brief
This file creates, set the behaviour and rules for the first boss Lancelot

Copyright (C) 2018 DigiPen Institute of Technology.
Reproduction or disclosure of this file or its contents
without the prior written consent of DigiPen Institute of
Technology is prohibited.
*/
/* End Header **************************************************************************/
#pragma once

#include "Characters.h"
#include "Dragon.h"

class Lancelot : public Characters
{
private:
	enum Lancelot_Moveset
	{
		// the attacks he will have
	};

	void Idle(float dt); // lancelot idle state

	void Moving(const Dragon &player, float dt); // moving state of lancelot

	void Attack(const Dragon &player, float dt); // attack state of lancelot

	void Init(void); //initializer for lancelot

	void Lancelot_Phase2(void);

	void (Lancelot::*ll_attacks[3])(void); // function pointer to lancelot attacks

	bool phase2;  // flag for lancelot second phase

public:
	Lancelot(); //default constructor for lancelot

	void Update() override {}; // ignored
	void Update(const Dragon &player, float dt);

};