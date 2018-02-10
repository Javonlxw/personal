/* Start Header ************************************************************************/
/*!
\file       Merlin.h
\author     Jacob Lim
\par email: jacob.lim\@digipen.edu
\brief
  Merlin class declared here.

Copyright (C) 20xx DigiPen Institute of Technology.
Reproduction or disclosure of this file or its contents
without the prior written consent of DigiPen Institute of
Technology is prohibited.
*/
/* End Header **************************************************************************/


#pragma once
#include "AEEngine.h"
#include "Characters.h"
#include "Boss_States.h"
#include "Dragon.h"
#include <vector>

//Attack state of Merlin
enum Merlin_Attack_State {
  MELEE = 0,    //Melee attack
  SINGLE_EBALL, //Single Shot Energy Ball
  SPREAD_EBALL, //Spread Shot Energy Ball
  ARROW_RAIN    //Arrow Rain
};

//Alias for Merlin Attack State
using MAS = Merlin_Attack_State;

const int    Merlin_HP       { 500 };     //Merlin's HP
const int    M_Phase2_HP     { 250 };     //Merlin's phase 2 HP
const int    Blink_CD_Time   { 360 };     //Cooldown time for Blink
const int    Merlin_Att_Inter{ 120 };     //Merlin's Attack Interval

const float  Eball_CD_Time   { 120.0f  }; //Energy Ball cooldown time
const float  Eball_Death     { 1000.0f }; //How far the energy ball travels

const float  Spread_CD_Time  { 240.0f  }; //Spread shot cooldown time
const float  Spread_Death    { 1000.0f }; //How far the spread shot travels

const int    A_Rain_Buffer   { 30      }; //How many arrows shot
const float  A_Rain_CD_Time  { 360.0f  }; //Arrow rain cooldown time
const float  A_Rain_Death    { 1000.0f }; //Arrow death

const float  Melee_CD_Time   { 60.0f };   //Melee cooldown time

const AEVec2 Blink_Pos_1     { -300.0f, 100.0f}; //(Point) 1st blink postion
const AEVec2 Blink_Pos_2     { 300.0f, 100.0f }; //(Point) 2nd blink postion
const AEVec2 Blink_Pos_3     { 0.0f, 180.0f   }; //(Point) 3rd blink postion

const float  Merlin_Start_X  { 320.0f  }; //Merlin default position X
const float  Merlin_Start_Y  { -120.0f }; //Merlin default position Y

class Merlin : public Characters {

public:
  Merlin();
  //Updates Merlin
  void Update(const Dragon &player);
  //Don't need to use dt
  void Update(float dt) { UNREFERENCED_PARAMETER(dt); }
  //Don't need default Update function
  void Update() {};
  //Renders Merlin and attacks
  void Render();
  //Returns boss state, use this to check if boss is attacking
  Boss_Action_State GetState() const { return M_Curr; }
  //Returns the type of attack
  MAS GetAttack() const { return M_Att_Curr; }
  //Returns Spread Eball to check for collision
  const std::vector<Boss_Attack> &GetSpreadBall() const { return Spread_Eball; }
  //Returns arrow to check for collision
  const std::vector<Boss_Attack> &GetArrow() const { return Arrow; }
  //Returns Eball to check for collision
  const Projectile &GetEball() const { return Eball; }
  ~Merlin();

private:
  void(Merlin::*Merlin_Attack)(const Dragon &player); //Pointer to current attack function
  void(Merlin::*Merlin_State)(const Dragon &player);  //Pointer to current state function
  void CheckState(const Dragon &player);              //Sets current state to next, checks for next state
  bool CheckAttack(const Dragon &player);             //Check if can attack, and sets the appropriate 
                                                      //function for it; Returns true if can attack, else false
  int castime;              //Arrow rain cast time
  int Attack_Interval;      //Time between attacks
  Boss_Action_State M_Curr; //Current state
  Boss_Action_State M_Next; //Next state
  MAS M_Att_Curr;           //Current attack state

  Boss_Attack M_Melee;                   //Melee attack
  Boss_Attack Eball;                     //Projectile for engery ball
  std::vector<Boss_Attack> Spread_Eball; //Projectiles for spread shot energy ball
  std::vector<Boss_Attack> Arrow;        //Projectiles for arrow rain
  //Blink attack struct
  struct Blink {
    int CD_Time  { Blink_CD_Time }; //Cooldown time
    bool Cooldown{     false     }; //Check if cooldown is active
  };
  Blink Blink_;       //Blink utilities
  Sprite MagicCircle; //Magic circle sprite
  Transform MC_Pos;   //Magic circle transform

  //Attack Functions START///////////////////////////////////////////////////////////////
  void Melee(const Dragon &player);    //Melee attack function
  void S_Eball(const Dragon &player);  //Single shot energy ball function
  void Sp_Eball(const Dragon &player); //Spread shot energy ball function
  void A_Rain(const Dragon &player);   //Arrow rain function
  //Attack Functions END/////////////////////////////////////////////////////////////////

  //State Functions START///////////////////////////////////////////////////////////////
  void Idle(const Dragon &player);   //Idle state function
  void Move(const Dragon &player);   //Move state function
  void Attack(const Dragon &player); //Attack state function
  //State Functions END/////////////////////////////////////////////////////////////////

};
