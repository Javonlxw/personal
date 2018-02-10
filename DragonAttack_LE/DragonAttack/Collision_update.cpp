///* Start Header ************************************************************************/
///*!
//\file       Collision_update.cpp
//\author     Javon Lee
//\par email: xiongweijavon.lee\@digipen.edu
//\brief
//This file contains the necessary functions to update the collision of each game object in each game state 
//
//Copyright (C) 2018 DigiPen Institute of Technology.
//Reproduction or disclosure of this file or its contents
//without the prior written consent of DigiPen Institute of
//Technology is prohibited.
//*/
///* End Header **************************************************************************/
//
//#include "Collision_update.h"
//
//void CollisionUpdate(GAME_STATE current)
//{
//	unsigned int i = 0, j = 0; //counter 
//	
//	//check if the current Game_state requires collision detection or not 
//	if (!GS_Collision(current))
//		return;
//	using namespace Col_update; 
//
//	for (i = 0; i < Max_Num_GameObj; ++i)
//	{
//		if (! GameObjList[i].IsActive())
//		{
//			continue; //skip to next iteration if the object is inactive 
//		}
//
//		if (GameObjList[i].Collision_.GetShape() == NoShape)
//		{
//			continue; //skip to next iteration if the object is a non-collidable object 
//		}
//
//		//Check if the Object is a Dragon 
//		if (GameObjList[i].Collision_.GetObjType() == Obj_Dragon)
//		{
//			//Check Dragon with other game objects 
//			for (j = 0; j < Max_Num_GameObj; ++j)
//			{
//				if (Check_Rectangle_Obj(GameObjList[j]))  //check if other object is a Rectangle 
//				{
//					if (GameObjList[i].Collision_.St_Rect_Rect(GameObjList[j].Collision_)) //perform static rectangle-rectangle collision check 
//					{
//						//There's COLLISION, input reaction. 
//					}
//
//					 //IF  Object's current MAGNITUDE/speed > 0, Currently idk where to access an object's magnitude 
//					if (GameObjList[i].Collision_.Dy_Rect_Rect(GameObjList[j].Collision_))
//					{
//						//There's COLLISION, input reaction. 
//					}
//
//				}// end of if-stmt for Rectangle 
//
//
//				if (Check_Circle_Obj(GameObjList[j]))
//				{
//					if (GameObjList[i].Collision_.St_Rect_Circle(GameObjList[j].Collision_)) //perform static rectangle-rectangle collision check 
//					{
//						//There's COLLISION, input reaction. 
//					}
//				}
//				
//
//
//
//			}
//		}
//
//	}
//
//}
//
//bool Check_Rectangle_Obj(GameObject& currObj)
//{
//	if (currObj.Collision_.GetShape == Rect)
//		return true;
//
//	return false;
//}
//
//bool Check_Circle_Obj(GameObject& currObj)
//{
//	if (currObj.Collision_.GetShape == Circle)
//		return true;
//
//	return false;
//}
//
//bool GS_Collision(GAME_STATE current)
//{
//	if (current == GS_CREDITS || current == GS_QUIT || current == GS_RESTART) // add in main menu after Engine Proof
//	{
//		return false;
//	}
//
//	return true;
//}