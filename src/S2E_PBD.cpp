
#ifdef WIN32
#define _WINSOCKAPI_    // stops windows.h including winsock.h
#include <tchar.h>
#include <windows.h>
#endif

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

// Simulator includes
#include "Interface/LogOutput/Logger.h"
#include "MCSimExecutor.h"

//Add custom include files
#include "Simulation/Case/PBD_Case.h"

// degub print of initialize file path
void print_path(std::string path)
{
#ifdef WIN32
  std::cout << path << std::endl;
#else
  const char *rpath = realpath(path.c_str(), NULL);
  if(rpath) {
    std::cout << rpath << std::endl;
    free((void *)rpath);
  }
#endif
}

// Main function
int main(int argc, char* argv[])
{
  //Set initialize file
  std::string ini_file = "../../data/ini/PBD_SimBase.ini";

  std::cout << "Starting simulation..." << std::endl;
  std::cout << "\tIni file: "; print_path(ini_file);

  auto simcase = PBD_Case(ini_file);
  simcase.Initialize();
  simcase.Main();

  return EXIT_SUCCESS;
}
