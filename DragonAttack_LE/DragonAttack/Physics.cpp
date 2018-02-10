#include "AEEngine.h"
#include "Physics.h"

float velocity = 5.0f;
float gravity = 0.5f;
float accel = 4.0f;

float ApplyMovement(float dt)
{
	return velocity + accel * dt;
}

float ApplyJump()
{
	return velocity * accel;
}

float ApplyGravity()
{
	return gravity;
}

float IfAirborne(float posY)
{
	if (posY <= 0)
	{
		accel = 4.0f;
	}
	else
	{
		accel -= gravity;
		posY += accel;
	}
	return posY;
}