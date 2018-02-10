#include "LE.h"
#include <iostream>

enum objType
{
  PLATFORM = 0,
  FLOOR,
  BG
};

AEGfxVertexList *mesh_plat;
AEGfxVertexList *mesh_floor;
AEGfxVertexList *mesh_bg;
AEGfxVertexList *grid;

int camX = 0, camY = 0;

struct oInst
{
  AEGfxVertexList* p_mesh = nullptr;
  objType type = PLATFORM;
  int x = 0;
  int y = 0;
  bool active = false;;
};

typedef struct TilePos {
  int x;
  int y;
}TP;

TP TileMap[16][16];

void UpdateTM()
{
  for (char i = 0; i < 16; ++i)
    for (char j = 0; j < 16; ++j)
    {
      TileMap[i][j].x = camX - 640 + (i * 128);
      TileMap[i][j].y = camY - 360 + (j * 72);
    }
}

void snap(int &curX, int &curY)
{
  int x = (curX / 128) * 128;
  int y = (curY / 72) * 72;
  for (char i = 0; i < 16; ++i)
    for (char j = 0; j < 16; ++j)
    {
      if (TileMap[i][j].x == x)
      {
        if (camX > 0)
          curX = TileMap[i][j].x + camX;
        else
          curX = TileMap[i][j].x - 640;
      }
      if (TileMap[i][j].y == y)
      {
        if(camY > 0)
          curY = TileMap[i][j].y - 360;
        else
          curY = 360 - TileMap[i][j].y;
      }
        
    }
}

std::vector <oInst> InstList;

oInst createInst(objType type, int x, int y)
{
  oInst tmp;
  switch (type)
  {
  case PLATFORM:
    tmp.p_mesh = mesh_plat;
    break;
  case FLOOR:
    tmp.p_mesh = mesh_floor;
    break;
  case BG:
    tmp.p_mesh = mesh_bg;
    break;
  default:
    tmp.p_mesh = nullptr;
    break;
  }
  tmp.type = type;
  tmp.x = x;
  tmp.y = y;
  tmp.active = true;
  return tmp;
}

namespace Level_Editor
{

  void Init(void)
  {
    mesh_plat = CreateRectangle(PLAT_WIDTH, PLAT_HEIGHT, 1.0f, 1.0f, 0xFF0000);
    mesh_floor = CreateRectangle(FLOOR_WIDTH, FLOOR_HEIGHT);
    mesh_bg = CreateRectangle(BG_WIDTH, BG_HEIGHT);
    grid = CreateSquare(1.0f, 1.0f, 1.0f, 0);
    UpdateTM();
  }

  void Load(void)
  {
    //Reserve memory space for 100 objects first
    InstList.reserve(100);
  }

  void Update(float dt)
  {
    UNREFERENCED_PARAMETER(dt);
    objType currType = PLATFORM;
    if (AEInputCheckTriggered(AEVK_P))
      currType = PLATFORM;
    if (AEInputCheckTriggered(AEVK_B))
      currType = BG;
    if (AEInputCheckTriggered(AEVK_F))
      currType = FLOOR;
    if (AEInputCheckTriggered(AEVK_LEFT))
      camX -= 128;
    if (AEInputCheckTriggered(AEVK_RIGHT))
      camX += 128;
    if (AEInputCheckTriggered(AEVK_UP))
      camY += 72;
    if (AEInputCheckTriggered(AEVK_DOWN))
      camY -= 72;
    //Update tile map
    UpdateTM();
    if (AEInputCheckTriggered(AEVK_LBUTTON))
    {
      int cursorX;
      int cursorY;
      AEInputGetCursorPosition(&cursorX, &cursorY);
      snap(cursorX, cursorY);
      std::cout << cursorX << " " << cursorY << std::endl;
      InstList.push_back(createInst(currType, cursorX, cursorY));
    }
  }

  void Draw(void)
  {
    AEGfxSetCamPosition((float)camX, (float)camY);
    for (unsigned i = 0; i < InstList.size(); ++i)
    {
        AEGfxSetRenderMode(AE_GFX_RM_COLOR);
        AEGfxSetPosition((float)InstList[i].x, (float)InstList[i].y);
        AEGfxTextureSet(nullptr, 0.0f, 0.0f);
        AEGfxSetTintColor(1.0f, 1.0f, 1.0f, 1.0f);
        AEGfxMeshDraw(InstList[i].p_mesh, AE_GFX_MDM_TRIANGLES);
    }
    for (char i = 0; i < 16; ++i)
      for (char j = 0; j < 16; ++j)
      {
        AEGfxSetRenderMode(AE_GFX_RM_COLOR);
        AEGfxSetPosition((float)TileMap[i][j].x, (float)TileMap[i][j].y);
        AEGfxTextureSet(nullptr, 0.0f, 0.0f);
        AEGfxSetTintColor(1.0f, 1.0f, 1.0f, 1.0f);
        AEGfxMeshDraw(grid, AE_GFX_MDM_TRIANGLES);        
      }
  }

  void Free(void)
  {

  }

  void Unload(void)
  {

  }

}