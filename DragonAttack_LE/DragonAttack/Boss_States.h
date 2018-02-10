#pragma once
#include <utility>
#include "Projectile.h"


// States for boss actions
enum Boss_Action_State
{
	IDLE = 0,
	MOVING,
	ATTACK,
	OBSTACLE
};



// A type to store boss attack and the cooldown
struct Boss_Attack : public Projectile
{
	float cooldown_timer;
	bool cooldown;

	Boss_Attack(Sprite&& t_attack, Col_Comp &&t_col)
		:cooldown_timer{0}, Projectile(std::move(t_attack), std::move(t_col)), 
		 cooldown{false}
	{
	}

	void Update() {}; // ignored
	void Update(float dt);  //has to check cooldown_timer
	void Pos(); 

};

//returns a random number within a range
int Get_Random_Num(int range);