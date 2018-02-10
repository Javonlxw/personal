/* Start Header ************************************************************************/
/*!
\file       Grunt.cpp
\author     William Yoong
\par email: william.yoong\@digipen.edu
\brief
Grunts to be created for the game.

Copyright (C) 20xx DigiPen Institute of Technology.
Reproduction or disclosure of this file or its contents
without the prior written consent of DigiPen Institute of
Technology is prohibited.
*/
/* End Header **************************************************************************/
#include "Grunt.h"
#include "Transform.h"
#include "Physics.h"
#include <iostream>

namespace // global variables just for THIS file
{
	const int grunt_hp = 30;
	bool PlayerSeen = false;
	bool PlayerInRange = false;
	float MovementX = 0.5f;
	float MovementY = 0.5f;
	float CurrXPos;
	float CurrYPos;
	float moveSpd = 3.0f;
}

Grunt::Grunt()
	: Characters(S_CreateSquare(100.0f, 1.0f, 1.0f, ".//Textures/grunt.png"),
		grunt_hp, { 0.0f, 0.0f, 5.0f, 5.0f, Rect })
{
	SetActive(true);
	PosX = 400.0f;
	PosY = -120.0f;
	CurrXPos = this->PosX;
	CurrYPos = this->PosY;
}

void Grunt::Update(float dt, const Dragon &d)
{
	UNREFERENCED_PARAMETER(dt);
	LineOfSight(d);

	if (PlayerSeen == true)
	{
		MoveTowardPlayer(d);

		if (PlayerInRange == true)
		{
			AttackPlayer(d);
		}
	}	
	else
	{
		Idle(d); // issue now is that if im returning back to idle after losing sight of player
				 // the currPos is still the same as the initial one set at the start
	}
}

void Grunt::Pos()
{
}

void Grunt::MoveTowardPlayer(const Dragon &d)
{
	/*
	s32 mx , my;
	AEInputGetCursorPosition(&mx, &my);
	mx -= 640; my -= 640;
	*/
		
	MovementX = (this->PosX) - d.PosX;
	MovementY = (this->PosY) - d.PosY;

	if (MovementX < 100.0f && MovementY < 50.0f &&
		MovementX > -100.0f && MovementY > -50.0f) // Attack Range (change accordingly)
	{
		PlayerInRange = true;
	}
	else
	{
		PlayerInRange = false;
	}
		
	PosX -= MovementX / 60;
	PosY -= MovementY / 60;
		
	this->Transform_.SetTranslate(PosX, PosY);
	this->Transform_.Concat();
}

void Grunt::LineOfSight(const Dragon &d)
{
	
	float playerDist = (d.PosX - this->PosX);
	if (playerDist <= 400.0f && playerDist >= -400.0f)  // vision range
	{
		PlayerSeen = true;
		std::cout << "Player Seen" << std::endl;
	}
	else
	{
		PlayerSeen = false;
	}
}

void Grunt::AttackPlayer(const Dragon &d)
{
	/*MovementX = (this->PosX) - d.PosX;
	MovementY = (this->PosY) - d.PosY;*/

	MovementX = d.PosX - (this->PosX);
	MovementY = d.PosY - (this->PosY);

	// if player moves out of range
	if (MovementX > 100.0f && MovementY > 50.0f ||
		MovementX < -100.0f && MovementY < -50.0f)
	{
		PlayerInRange = false;
	}
	else
	{
		PlayerInRange = true;
	}
	//attack animation
	std::cout << "Attacking" << std::endl;
}

void Grunt::Idle(const Dragon &d)
{
	UNREFERENCED_PARAMETER(d);

	float MaxXPos = CurrXPos + 100.0f; // max boundary
	float MinXPos = CurrXPos - 100.0f; // min boundary

	if (PlayerSeen == false)
	{
		if ((PosX >= MaxXPos) || (PosX <= MinXPos))
		{
			moveSpd *= -1;
		}

		PosX += moveSpd;
		this->Transform_.SetTranslate(PosX, PosY);
		this->Transform_.Concat();
	}
}