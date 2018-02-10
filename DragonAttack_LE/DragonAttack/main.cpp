#include "AEEngine.h"
#include "Sprite.h"
#include "Transform.h"
#include "Create_Object.h"
#include "GameObject.h"
#include "Camera.h"
#include "GameStateManager.h"
#include "Collision.h"
#include <utility>
#include <iostream>

int WINAPI WinMain(HINSTANCE instanceH, HINSTANCE prevInstanceH, LPSTR command_line, int show)
{
  UNREFERENCED_PARAMETER(instanceH);
  UNREFERENCED_PARAMETER(prevInstanceH);
  UNREFERENCED_PARAMETER(command_line);
  UNREFERENCED_PARAMETER(show);

  const int R_WIDTH{ 1280 }; //Screen width
  const int R_HEIGHT{ 720 }; //Screen height

  float dt = 0.0f;       // delta time
  float timePast = 0.0f; // time past in game

  AESysInit(instanceH, show, R_WIDTH, R_HEIGHT, true, 60, NULL);
  AESysSetWindowTitle("My New Demo!");
  AESysReset();
  GSM::GSM_Init();

  while (GSM::IsGameRunning())
 {
	GSM::Init_and_Load();

    while (GSM::current == GSM::next)
    {
      AESysFrameStart();
      GSM::Update_and_Draw(dt);
      AESysFrameEnd();

	  dt = (f32)AEFrameRateControllerGetFrameTime();
	  timePast += dt;
    }



    GSM::Cleanup();
  }
  AESysExit();
  return 0;
}
