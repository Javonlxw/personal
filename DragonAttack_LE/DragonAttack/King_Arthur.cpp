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
#include <iostream>
#include "King_Arthur.h"
#include "Boss_States.h"
#include <vector>

namespace {

	std::vector <Boss_Attack> arthur; //an array of boss attacks

	std::vector <Characters> spawn_mobs; //an array to store the mobs to be spawn

	Boss_Action_State current_action = ATTACK; // different states of boss arthur

	const int health = 3000; // initial hp for king arthur

	const int phase2_hp = 1500; //phase 2 trigger

	const int interval = 50; // interval for triple slash

	const int range_limit = 300; // range limit for slash

	const char limit = 5; // num of king arthur attacks

	char behavior_swap = 0; // switch between idling and moving

	void Move_KA(float dt, King_Arthur &ka);

	void Set_Attack_Dir(King_Arthur &ka);
	
}



King_Arthur::King_Arthur(void)
	: Characters(S_CreateSquare(100.0f, 1.0f, 1.0f, ".//Textures/download.png"),
	  health, { 0.0f, 0.0f, 5.0f, 5.0f, Rect } ), phase2{false}
{
	(void) this->Transform_.SetTranslate(200.0f, Start_Pos_Y);
	this->Transform_.Concat(); // spawn king arthur at the location set
	this->PosX =  200.0f;      // change king arthur coordinates to the location set
	this->PosY = Start_Pos_Y;  // change king arthur coordinates to the location set
	this->SetActive(true);     // show him on screen
	this->Set_Direction(LEFT); // set king arthur to face left at the start
	this->SetVelocity({ 120, 0 });    // velocity for king arthur
	this->Sprite_.SetAlphaTransBM(1.0f, 1.0f, AE_GFX_BM_BLEND);
	Init();  // call initializer for king arthur move set
}

void King_Arthur::Init(void)
{
	const char* tex_ss = ".//Textures/1.jpg"; 
	const char* tex_3s = ".//Textures/holylight.jpg";
	const char* tex_ja = ".//Textures/light_texture2320.jpg";

	arthur.reserve(3);

	arthur.emplace_back(Boss_Attack(S_CreateSquare(20.0f, 1.0f, 1.0f, tex_ss), // for single slash
					 Col_Comp(0.0f, 0.0f, 5.0f, 5.0f, Rect))); 

	arthur.emplace_back(Boss_Attack(S_CreateSquare(10.0f, 1.0f, 1.0f, tex_3s), // 1st triple slash 
					 Col_Comp(0.0f, 0.0f, 5.0f, 5.0f, Rect)));

	arthur.emplace_back(Boss_Attack(S_CreateSquare(10.0f, 1.0f, 1.0f, tex_3s), // 2nd single slash
		Col_Comp(0.0f, 0.0f, 5.0f, 5.0f, Rect)));

	arthur.emplace_back(Boss_Attack(S_CreateSquare(10.0f, 1.0f, 1.0f, tex_3s), // 3rd triple slash
		Col_Comp(0.0f, 0.0f, 5.0f, 5.0f, Rect)));

	arthur.emplace_back(Boss_Attack(S_CreateSquare(30.0f, 1.0f, 1.0f, tex_ja), // jump attack
					 Col_Comp(0.0f, 0.0f, 5.0f, 5.0f, Rect)));

	(void)arthur[0].SetVelocity(AEVec2{ 3.0f, 0.0f }); // velocity of the slashes
	(void)arthur[1].SetVelocity(AEVec2{ 3.0f, 0.0f }); // velocity of the slashes
	(void)arthur[2].SetVelocity(AEVec2{ 3.0f, 0.0f }); // velocity of the slashes
	(void)arthur[3].SetVelocity(AEVec2{ 3.0f, 0.0f }); // velocity of the slashes
														 
	arthur[0].Transform_.SetScale(2.0f, 2.0f); // determine the size of projectile
	arthur[0].Transform_.Concat();

	for (char i = 0; i < limit; ++i) // spawning location of slashes
		arthur[i].PosY = -120.0f;

	ka_attacks[SINGLE_SLASH]    = &King_Arthur::Single_Slash;
	ka_attacks[TRIPLE_SLASH]    = &King_Arthur::Triple_Slash;
	ka_attacks[UNIQUE_MECHANIC] = &King_Arthur::Jump_Attack;

}

void King_Arthur::Update(float dt, const Dragon &d) 
{

	if (this->Get_HP() < phase2_hp) // activate phase 2 once hp drops is 50%
	{
		King_Arthur_Phase2();
	}

	if (behavior_swap == 3)
	{
		current_action = MOVING;  // switch state to moving
		behavior_swap = 0;        // reset the counter
	}

	// switch between the boss states
	switch (current_action)
	{
	case IDLE: this->Idle(dt);
		break;

	case MOVING: this->Moving(d,dt);
		break;

	case OBSTACLE: this->AvoidingObstacle();
		break;

	case ATTACK: this->Attack(d, dt);
		break;
	}

	for (int i = 0; i < limit; ++i) // update cooldowns on attacks
			arthur[i].Update(dt);


	
}

/******************************************************************************/
/*!
\fn         Idle

\brief      King arthur pauses for awhile for the player to refocus

\param      float dt
- delta time to calculate how long king arthur freezes

\retval     void
No return.
*/
/******************************************************************************/
void King_Arthur::Idle(float dt)
{
	UNREFERENCED_PARAMETER(dt);
	static float idle = 2.0f; //determine how long king arthurs stops

	if (!idle) //sets the timer to 3 seconds 
	{
		idle = 2.0f;
	}
	
	//changes state to attack once idling is finished, reset idle
	idle <= 0? current_action = ATTACK, idle = 0 : idle -= dt;
	Set_Attack_Dir(*this);
	
}

void King_Arthur::Moving(const Dragon &d, float dt)
{
	UNREFERENCED_PARAMETER(dt); // need to check for platform in phsae 2
	
	static float move_duration = 3.0f; // duration KA moves
	
	if (!move_duration) //set the timer for moving
		move_duration = 3.0f;
	
	if (d.PosX - this->PosX == 0) // character and boss intersect
		this->Set_Direction(STAY);
	else
	//determine if dragon is right or left of king arthur
	(d.PosX - this->PosX) > 0 ? this->Set_Direction(RIGHT):
	this->Set_Direction(LEFT);

	/*if(phase2)
		d.PosY - this->PosY;*/
	
	Move_KA(dt, *this);

	//change state to attack once move_duration is 0, reset move_duration
	move_duration <= 0 ? current_action = ATTACK, move_duration = 0 :
		                 move_duration -= dt;

	Set_Attack_Dir(*this);
}

void King_Arthur::Attack(const Dragon &d, float dt)
{
	UNREFERENCED_PARAMETER(dt);
	UNREFERENCED_PARAMETER(d);

	KA_MoveSet currAttk = SINGLE_SLASH; //default attack for king arthur
	
	if (d.PosX - this->PosX == 0) // character and boss intersect
		this->Set_Direction(STAY);
	else
		//determine if dragon is right or left of king arthur
		(d.PosX - this->PosX) > 0 ? this->Set_Direction(RIGHT) :
		this->Set_Direction(LEFT);

	//unique mechanic has the highest priority 
	if (! (arthur[4].cooldown) )
		currAttk = UNIQUE_MECHANIC;

	//followed by triple slash
	else if (! (arthur[1].cooldown) ) 
		currAttk = TRIPLE_SLASH; 

	(this->*ka_attacks[currAttk])(); // call the coresponding attack function
}

void King_Arthur::AvoidingObstacle(void)
{
	//move to player in phase 2 with platforms
	//path finding i guess?
}

void King_Arthur::King_Arthur_Phase2(void)
{
	// change the phase 1 machanic to phase 2
	// play some cinematic shit
	// healing should be just particle effects

	arthur.pop_back(); // delete jump attack

	const char* tex_ja = ".//Textures/light_texture2320.jpg";
	arthur.emplace_back(Boss_Attack(S_CreateSquare(30.0f, 1.0f, 1.0f, tex_ja),
		Col_Comp(0.0f, 0.0f, 5.0f, 5.0f, Rect))); // add the 2nd mechanic, still in progress

	this->phase2 = true; // set phase 2 flag to true

	// change unique mechanism pointer to heal and spawn for phase 2
	ka_attacks[UNIQUE_MECHANIC] = &King_Arthur::Heal_and_Spawn; 
}

void King_Arthur::Jump_Attack()
{
	if (arthur[4].cooldown) // skill still on cooldown
		return;

	static const int left_boundary = -600; // boundaries of charge attack
	static const int right_boundary = 600; // boundaries of charge attack

	this->SetVelocity({ 600, 0 }); // change king arthur's velocity to dash
	Move_KA(0.016f, *this); // i am sorry for this

	// need check for collision against dragon
	if(this->PosX > right_boundary || this->PosX < left_boundary)
	{
		current_action = IDLE;           // switch current state to idle
		arthur[4].cooldown = true;       // start the cooldown of dash attk
		this->SetVelocity({ 120, 0 });   // reset the velocity of king arthur
		++behavior_swap;                 
	}

	this->Transform_.Concat();

	arthur[4].cooldown_timer = 10.0f; // cooldown of jump attacl
}

void King_Arthur::Single_Slash(void)
{
	if (!(arthur[0].cooldown_timer)) // sets the attack to start from KA
	{
		arthur[0].PosX = this->PosX;
		arthur[0].SetActive(true);
	}
	
	if (arthur[0].GetDist() > range_limit) //disappears once it flys 500 unit away
	{
		arthur[0].SetActive(false); // makes it disappears form screen

		arthur[0].ResetDist(); // error here, need to reset it to king arthur's position			

		current_action = IDLE; // set the behaviour to idle

		arthur[0].PosX = 0.0f;  // to remove flicker

		arthur[0].cooldown = true; // skill on cooldown

		++behavior_swap;

	}

	arthur[0].Projectile::Pos(); // need to check what the rotation does
	arthur[0].Projectile::Update(); // check what the sqrt does also
	arthur[0].Render();

	arthur[0].cooldown_timer = 2.0f;
}

void King_Arthur::Triple_Slash(void)
{
	char ts_limit = limit - 1; // limit the loop to just the triple slashes

	if (arthur[1].cooldown) // skill still on cooldown
		return;
	
	if (!(arthur[1].cooldown_timer)) // sets the attack to start from KA
	{
		for (char i = 1; i < ts_limit; ++i)
		{
			arthur[i].SetActive(true);
			arthur[i].PosX = this->PosX;
		}
	}

	//slashes 3 times at the player with an interval inbetween
	for (char i = 1; i < ts_limit; ++i)
	{
		if (arthur[i].GetDist() > range_limit) // makes projectiles disappear
		{
			arthur[i].SetActive(false); // disappears from screen
			arthur[i].ResetDist();     // reset the projectile distance traveled
			arthur[i].PosX = 0.0f; // remove flicker

			if (i == 3) // last slash has completed
			{
				current_action = IDLE; // change state to idle
				arthur[1].cooldown = true; // start cooldown		
				++behavior_swap;
			}
		}
	}
	// issue here if i want to use jacob projectile update and pos
	if (arthur[1].GetDist() > interval || !(arthur[1].IsActive()))
	{
		//fire the second projectile once the first 1 travels 50 pixels
		arthur[2].Pos();

		if (arthur[2].GetDist() > interval || !(arthur[2].IsActive()))
			// fires the third one once the second 1 travels 0 pixel
			arthur[3].Pos();
	}

	arthur[1].Pos(); // render the first slash
	arthur[1].cooldown_timer = 10.0f;


}

void King_Arthur::Heal_and_Spawn(void)
{
	//teleport to a platform, heal and spawn mobs

	if (!(arthur[4].cooldown_timer)) // teleport once
	{
		//get a random number and teleport to the coresponding platform
		switch (Get_Random_Num(3)) 
		{
		case 1: this->PosX = 300.0f, this->PosY = 100.0f;
			break;

		case 2: this->PosX = -300.0f, this->PosY = 100.0f;
			break;

		case 3: this->PosX = 50.0f, this->PosY = 200.0f;
			break;

		default: break;
		}
		this->Transform_.SetTranslate(this->PosX, this->PosY);
		this->Transform_.Concat();
	}

	arthur[4].cooldown_timer = 20.0f;
}

King_Arthur::~King_Arthur(void)
{
	arthur.clear();
}


namespace
{
	void Move_KA(float dt, King_Arthur &ka)
	{
		if (ka.Get_Direction() == RIGHT) // set all attacks to go right
		{
			ka.PosX += ka.GetVelocity().x * dt; // move ka to the right
			(void)ka.Transform_.SetTranslate(ka.PosX, ka.PosY);
		}
		else if (ka.Get_Direction() == LEFT) // set all attacks to go left
		{
			ka.PosX -= ka.GetVelocity().x * dt; // move KA to the left
			(void)ka.Transform_.SetTranslate(ka.PosX, ka.PosY);
		}

		ka.Transform_.Concat();
	}

	void Set_Attack_Dir(King_Arthur & ka)
	{
		
		if (ka.Get_Direction() == RIGHT) // set all attacks to go right
		{
			for (int i = 0; i < limit; ++i)
			{
				arthur[i].SetDir(true);
			}
		}
		else if (ka.Get_Direction() == LEFT) // set all attacks to go left
		{
			for (int i = 0; i < limit; ++i)
			{
				arthur[i].SetDir(false);
			}
		}
	}
}

