#include "OBC_Sat1.h"
#include "../../Simulation/Spacecraft/PBD_Components.h"

OBC_Sat1::OBC_Sat1(ClockGenerator* clock_gen, PBD_Components& components) : OBC(clock_gen), components_(components)
{
  Initialize();
}

OBC_Sat1::~OBC_Sat1()
{
}

void OBC_Sat1::Initialize()
{
}

void OBC_Sat1::MainRoutine(int count)
{
}
