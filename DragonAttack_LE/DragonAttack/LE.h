#pragma once
#include "AEEngine.h"
#include "GameObject.h"
#include "Create_Object.h"
#include <cmath>
#include <vector>

namespace Level_Editor
{

  void Init(void);

  void Load(void);

  void Update(float dt);

  void Draw(void);

  void Free(void);

  void Unload(void);

}