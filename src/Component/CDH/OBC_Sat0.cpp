#include "OBC_Sat0.h"
#include "../../Simulation/Spacecraft/PBD_Components.h"

OBC_Sat0::OBC_Sat0(ClockGenerator* clock_gen, PBD_Components& components) : OBC(clock_gen), components_(components)
{
  Initialize();
}

OBC_Sat0::~OBC_Sat0()
{
}

void OBC_Sat0::Initialize()
{
}

void OBC_Sat0::MainRoutine(int count)
{
}
