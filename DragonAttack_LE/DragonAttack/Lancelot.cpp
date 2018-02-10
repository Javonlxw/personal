#include "Lancelot.h"
#include "Boss_States.h"
#include <iostream>
#include <vector>

/*
mad enhancement 
move faster
faster attack speed 

phase 2
arondight longer reach?
*/


namespace
{
	std::vector<Boss_Attack> lancelot; // array to store lancelot attack

	Boss_Action_State current_action = MOVING;

	const int health = 1500;

}

Lancelot::Lancelot(void)
	: Characters(S_CreateSquare(100.0f, 1.0f, 1.0f, ".//Textures/lancelot.jpg"),
	  health, { 0.0f, 0.0f, 5.0f, 5.0f, Rect }), phase2{ false } 
{
	(void) this->Transform_.SetTranslate(200.0f, Start_Pos_Y); // spawn lancelot at this location
	this->Transform_.Concat();
	this->PosX = 200.0f;      // update lancelot current coordinates
	this->PosY = Start_Pos_Y; // update lancelot current coordinates
	this->SetActive(true);    // spawn lancelot
	this->Set_Direction(LEFT); 
	this->Sprite_.SetAlphaTransBM(1.0f, 1.0f, AE_GFX_BM_BLEND);

}

void Lancelot::Init()
{
	//initialize lancelot attacks here
}

void Lancelot::Idle(float dt)
{
	UNREFERENCED_PARAMETER(dt);
}

void Lancelot::Attack(const Dragon &d, float dt)
{
	UNREFERENCED_PARAMETER(d);
	UNREFERENCED_PARAMETER(dt);
}

void Lancelot::Moving(const Dragon &d, float dt)
{
	UNREFERENCED_PARAMETER(d);
	UNREFERENCED_PARAMETER(dt);
}

void Lancelot::Update(const Dragon &d, float dt)
{
	UNREFERENCED_PARAMETER(d);
	UNREFERENCED_PARAMETER(dt);
}

void Lancelot::Lancelot_Phase2(void)
{
	this->phase2 = true;
}