/* Start Header ************************************************************************/
/*!
\file       Input_Handler.h
\author     William Yoong
\par email: william.yoong\@digipen.edu
\brief
Input from user handled here.

Copyright (C) 20xx DigiPen Institute of Technology.
Reproduction or disclosure of this file or its contents
without the prior written consent of DigiPen Institute of
Technology is prohibited.
*/
/* End Header **************************************************************************/

#pragma once
#include "GameStateManager.h"
#include "AEEngine.h"
#include <vector>

//variables to store the different keys being use in the game
namespace Input
{
	const int jump_up    = AEVK_SPACE;
	const int go_left    = AEVK_A;
	const int go_right   = AEVK_D;
	const int go_up      = AEVK_W;
	const int go_down    = AEVK_S;
	const int fire       = AEVK_RETURN;
	const int special    = AEVK_P;
	const int quit_game  = AEVK_ESCAPE;


	/******************************************************************************/
	/*!
	\fn         void Get_User_Input

	\brief      Get the key(s) that the user pressed

	\param      void

	\retval     void
	No return.
	*/
	/******************************************************************************/
	const std::vector <int> &Get_User_Input(void);

	void ClearBuffer(void);

}