/* Start Header ************************************************************************/
/*!
\file       Grunt.h
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

#pragma once
#include "Characters.h"
#include "Dragon.h"

class Grunt : public  Characters
{
private: 
	void LineOfSight(const Dragon &d);

	void MoveTowardPlayer(const Dragon &d);

	void AttackPlayer(const Dragon &d);

	void Idle(const Dragon &d);
public:
	/******************************************************************************/
	/*!
	\fn         Grunt

	\brief      Constructor for a Grunt Object

	\param      Sprite&& t_sprite
					- an rvalue of the sprite to be moved

	\retval     void
	No return.
	*/
	/******************************************************************************/
	Grunt();

	/******************************************************************************/
	/*!
	\fn         Grunt

	\brief      Constructor for a Grunt Object

	\param      Sprite&& t_sprite
					- an rvalue of the sprite to be moved
				x 
					- location to spawn the enemy
				y 
					- location to spawn the enemy

	\retval     void
	No return.
	*/
	/******************************************************************************/
	//Grunt(Sprite&& t_sprite, const float x, const float y);

	void Update(float dt, const Dragon &d);
	void Update() override {}; // ignored

	void Pos();

	//~Grunt(); // destructor
};

/******************************************************************************/
/*!
\fn         void *Create_Grunt

\brief      Creates a grunt enemy

\param      void

\retval     void
No return.
*/
/******************************************************************************/
//Grunt *Create_Grunt();

/******************************************************************************/
/*!
\fn         void *Delete_Grunt

\brief      Deletes a grunt enemy

\param      A pointer to a grunt object

\retval     void
No return.
*/
/******************************************************************************/
//void Delete_Grunt(Grunt* grunt);
