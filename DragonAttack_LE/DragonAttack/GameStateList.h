/* Start Header ************************************************************************/
/*!
\file GameStateList.h
\author William Yoong, william.yoong, 1701674

\par william.yoong\@digipen.edu)
\date Jan 24, 2017
\brief
This header file contains the enum for the different states the game haves

Copyright (C) 20xx DigiPen Institute of Technology.
Reproduction or disclosure of this file or its contents
without the prior written consent of DigiPen Institute of
Technology is prohibited.
*/
/* End Header **************************************************************************/
#pragma once

/* The different game states for the game */
enum GAME_STATE
{
  LE = 0,
	GS_MAIN,
	GS_LEVEL1,
	GS_LEVEL2,
	GS_LEVEL3,

	GS_CREDITS,
	GS_QUIT,
	GS_RESTART
};