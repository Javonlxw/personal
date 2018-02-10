#include "Boss_States.h"
#include <cstdlib>
#include <ctime>
#include <iostream>

void Boss_Attack::Update(float dt)  //has to check cooldown_timer
{
	if (cooldown)
	{
		cooldown_timer -= dt;

		if (cooldown_timer <= 0)
		{
			cooldown_timer = 0;
			cooldown = false;
		}
	}


}
void Boss_Attack::Pos()
{
	//calculate the position and shit for attack
	
	if (! IsActive()) // only render when its shown on screen
		return;
	
	if (this->GetDir()) // shoots projectile to the right
		this->PosX += this->GetVelocity().x;
	else // shoots projectile to the left
		this->PosX -= this->GetVelocity().x;

	//moves the projectile in the correct direction frame by frame
	this->AddDist (this->GetVelocity().x);
	(void)this->Transform_.SetScale(2.0f,2.0f);
	(void)this->Transform_.SetTranslate(this->PosX, this->PosY);
	this->Transform_.Concat();
	
	// show it on screen
	this->Render();

}

int Get_Random_Num(int range)
{
	srand((unsigned) time(NULL));

	return rand() % range + 1;
}