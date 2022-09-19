## General Settings
cmake_policy(SET CMP0048 NEW)
# set(CMAKE_CXX_STANDARD 14)
project(S2E_PBD
    LANGUAGES CXX
    DESCRIPTION "S2E_PBD"
    VERSION 1.0
)

cmake_minimum_required(VERSION 3.13)

# preprocessor
if(WIN32)
	add_definitions(-DWIN32)
endif()

## set directory path
set(S2E_CORE_DIR ../s2e-core_pbd)
if(NOT DEFINED EXT_LIB_DIR)
  set(EXT_LIB_DIR ../ExtLibraries)
endif()
set(CSPICE_DIR ${EXT_LIB_DIR}/cspice)
set(NRLMSISE00_DIR ${EXT_LIB_DIR}/nrlmsise00)
set(EIGEN_DIR ${EXT_LIB_DIR}/Eigen)

if(NOT DEFINED FLIGHT_SW_DIR)
  set(FLIGHT_SW_DIR ../FlightSW)
endif()
if(NOT DEFINED C2A_NAME)
  set(C2A_NAME "c2a_oss")
endif()
set(FLIGHT_SW_DIR ../FlightSW)

if(WIN32)
elseif(APPLE)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ./out)  #set build path
endif()

## options to use C2A
option(USE_C2A "Use C2A")
set(USE_C2A OFF)
if(USE_C2A)
  add_definitions(-DUSE_C2A)
  add_definitions(-DSILS_FW)
  #include_directories
  include_directories(${C2A_DIR}/SH7254R_C2A/src_core/System/AnomalyLogger)
  include_directories(${C2A_DIR}/SH7254R_C2A/src_core/System/ApplicationManager)
  include_directories(${C2A_DIR}/SH7254R_C2A/src_core/System/ModeManager)
  include_directories(${C2A_DIR}/SH7254R_C2A/src_core/System/TaskManager)
  include_directories(${C2A_DIR}/SH7254R_C2A/src_core/System/WatchdogTimer)
  include_directories(${C2A_DIR}/SH7254R_C2A/src_core/System/TimeManager)
  include_directories(${C2A_DIR}/SH7254R_C2A/src_core/CmdTlm)
  include_directories(${C2A_DIR}/SH7254R_C2A/src_user/Applications)
  include_directories(${S2E_CORE_DIR}/src/Interface/SpacecraftInOut)
  #add subdirectory
  add_subdirectory(${C2A_DIR} C2A_CORE)
  message("USE C2A")
endif()

set(S2E_DIR ${CMAKE_CURRENT_COURCE_DIR})

## include directories of S2E-CORE
include_directories(${CSPICE_DIR}/include)
include_directories(${NRLMSISE00_DIR}/src)
include_directories(${S2E_CORE_DIR}/src)
include_directories(${S2E_CORE_DIR}/src/Library/math)
include_directories(${S2E_CORE_DIR}/src/Library/utils)
include_directories(${S2E_CORE_DIR}/src/Dynamics)
include_directories(${S2E_CORE_DIR}/src/Disturbance)
include_directories(${S2E_CORE_DIR}/src/Interface/InitInput)
include_directories(${S2E_CORE_DIR}/src/Interface/LogOutput)
include_directories(${S2E_CORE_DIR}/src/Interface/SpacecraftInOut/Ports)
include_directories(${S2E_CORE_DIR}/src/Simulation)
include_directories(${S2E_CORE_DIR}/src/Simulation/MCSim)
include_directories(${S2E_CORE_DIR}/src/Simulation/Spacecraft)
include_directories(${S2E_CORE_DIR}/src/Simulation/Spacecraft/Structure)
include_directories(${S2E_CORE_DIR}/src/Simulation/GroundStation)
include_directories(${S2E_CORE_DIR}/src/Simulation/Case)
include_directories(${S2E_CORE_DIR}/src/Disturbance)
include_directories(${S2E_CORE_DIR}/src/Environment/Global)
include_directories(${S2E_CORE_DIR}/src/Environment/Local)
include_directories(${S2E_CORE_DIR}/src/Component)
include_directories(${S2E_CORE_DIR}/src/Component/Abstract)
include_directories(${S2E_CORE_DIR}/src/Component/AOCS)
include_directories(${S2E_CORE_DIR}/src/Component/CDH)
include_directories(${S2E_CORE_DIR}/src/Component/CommGS)
include_directories(${S2E_CORE_DIR}/src/Component/Mission/Telescope)
include_directories(${S2E_CORE_DIR}/src/Component/Power)
include_directories(${S2E_CORE_DIR}/src/Component/Propulsion)
include_directories(${S2E_CORE_DIR}/src/Component/Thermal)
include_directories(${EIGEN_DIR})
include_directories(${S2E_CORE_DIR}/src/Simulation/InterSatComm)
include_directories(${S2E_CORE_DIR}/src/Simulation/InterSatComm/Abstract)
include_directories(${S2E_CORE_DIR}/src/RelativeInformation)

## add_subdirectories
add_subdirectory(${S2E_CORE_DIR}/src/Component S2E_CORE/Component)
add_subdirectory(${S2E_CORE_DIR}/src/Disturbance S2E_CORE/Disturbance)
add_subdirectory(${S2E_CORE_DIR}/src/Dynamics S2E_CORE/Dynamics)
add_subdirectory(${S2E_CORE_DIR}/src/Environment/Global S2E_CORE/Environment/Global)
add_subdirectory(${S2E_CORE_DIR}/src/Environment/Local S2E_CORE/Environment/Local)
add_subdirectory(${S2E_CORE_DIR}/src/RelativeInformation S2E_CORE/RelativeInformation)
add_subdirectory(${S2E_CORE_DIR}/src/Interface/InitInput S2E_CORE/InitInput)
add_subdirectory(${S2E_CORE_DIR}/src/Interface/LogOutput S2E_CORE/LogOutput)
add_subdirectory(${S2E_CORE_DIR}/src/Interface/SpacecraftInOut S2E_CORE/SpacecraftInOut)
add_subdirectory(${S2E_CORE_DIR}/src/Interface/HilsInOut S2E_CORE/HilsInOut)
add_subdirectory(${S2E_CORE_DIR}/src/Library/igrf S2E_CORE/igrf)
add_subdirectory(${S2E_CORE_DIR}/src/Library/inih S2E_CORE/inih)
add_subdirectory(${S2E_CORE_DIR}/src/Library/math S2E_CORE/math)
add_subdirectory(${S2E_CORE_DIR}/src/Library/nrlmsise00 S2E_CORE/nrlmsise00)
add_subdirectory(${S2E_CORE_DIR}/src/Library/sgp4 S2E_CORE/sgp4)
add_subdirectory(${S2E_CORE_DIR}/src/Library/utils S2E_CORE/utils)
add_subdirectory(${S2E_CORE_DIR}/src/Library/optics S2E_CORE/optics)
add_subdirectory(${S2E_CORE_DIR}/src/Library/RelativeOrbit S2E_CORE/RelativeOrbit)
add_subdirectory(${S2E_CORE_DIR}/src/Library/Orbit S2E_CORE/Orbit)
add_subdirectory(${S2E_CORE_DIR}/src/Library/Geodesy S2E_CORE/Geodesy)
add_subdirectory(${S2E_CORE_DIR}/src/Simulation S2E_CORE/Simulation)

set(SOURCE_FILES
  src/S2E_PBD.cpp
  src/Simulation/Case/PBD_Case.cpp
  src/Simulation/Spacecraft/PBD_Sat.cpp
  src/Simulation/Spacecraft/PBD_Components.cpp
  src/Simulation/GroundStation/PBD_GroundStation.cpp
  src/Simulation/GroundStation/PBD_GSComponents.cpp
  src/Simulation/InterSatComm/PBD_InterSatComm.cpp
  src/Simulation/InterSatComm/RFSystem/RFSystemBeam.cpp
  src/Simulation/PBD/PBD_dgps.cpp
  src/Simulation/PBD/PBD_GnssObservation.cpp
  src/Simulation/PBD/PBD_Lambda.cpp
  src/Simulation/PBD/PBD_const.cpp
  src/Component/CDH/OBC_Sat0.cpp
  src/Component/CDH/OBC_Sat1.cpp
  src/Component/AOCS/PBD_GNSSReceiver.cpp
  src/Component/AOCS/InitGnssReceiver.cpp
  src/Component/RFSystem/RFSystemReceiver.cpp
  src/Component/RFSystem/RFSystemTransmitter.cpp
  src/Interface/InitInput/InitComponent/InitRFSystemReceiver.cpp
  src/Interface/InitInput/InitComponent/InitRFSystemTransmitter.cpp
  src/Interface/InitInput/InitInterSatComm/InitRFSystemBeam.cpp
)
# Create executable file
add_executable(${PROJECT_NAME} ${SOURCE_FILES})
# Compile option
if(MSVC)
  target_compile_options(${PROJECT_NAME} PUBLIC "/W4")
  target_compile_options(${PROJECT_NAME} PUBLIC "/source-charset:utf-8")
else()
  target_compile_options(${PROJECT_NAME} PUBLIC "${CMAKE_CXX_FLAGS}-Wall")
  set(CMAKE_CXX_FLAGS "-m32 -rdynamic -Wall -g")
  set(CMAKE_C_FLAGS "-m32 -rdynamic -Wall -g")
endif()

## cspice library
if(CYGWIN)
  SET (CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")
  find_library(CSPICE_LIB
    NAMES cspice.a csupport.a
    PATHS ${CSPICE_DIR}/cspice_cygwin/lib)
elseif(UNIX)
  find_library(CSPICE_LIB
    NAMES cspice.a csupport.a
    PATHS ${CSPICE_DIR}/cspice_unix/lib)
elseif(WIN32)
  SET (CMAKE_FIND_LIBRARY_SUFFIXES ".lib")
  find_library(CSPICE_LIB
    NAMES cspice.lib csupport.lib
    PATHS ${CSPICE_DIR}/cspice_msvs/lib)
endif()
if(NOT CSPICE_LIB)
  message(FATAL_ERROR "cspice not found in ${EXT_LIB_DIR}")
endif()

## nrlmsise00 library
if(CYGWIN)
  SET (CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  find_library(NRLMSISE00_LIB 
    NAMES libnrlmsise00.a
    PATHS ${NRLMSISE00_DIR}/lib)
elseif(UNIX)
  find_library(NRLMSISE00_LIB 
    NAMES libnrlmsise00.a
    PATHS ${NRLMSISE00_DIR}/lib)
elseif(WIN32)
  SET (CMAKE_FIND_LIBRARY_SUFFIXES ".lib")
  find_library(NRLMSISE00_LIB 
    NAMES libnrlmsise00.lib
    PATHS ${NRLMSISE00_DIR}/lib)
endif()
if(NOT NRLMSISE00_LIB)
  message(FATAL_ERROR "nrlmsise00 not found in ${EXT_LIB_DIR}")
endif()

## Linking libraries
set(S2E_LIBRARIES
  IGRF WRAPPER_NRLMSISE00 INIH SGP4 UTIL OPTICS RELATIVE_ORBIT_MODELS ORBIT_MODELS GEODESY MATH
)

# Initialize link
target_link_libraries(COMPONENT DYNAMICS GLOBAL_ENVIRONMENT LOCAL_ENVIRONMENT SC_IO RELATIVE_INFO ${S2E_LIBRARIES})
target_link_libraries(DYNAMICS GLOBAL_ENVIRONMENT LOCAL_ENVIRONMENT SIMULATION ${S2E_LIBRARIES})
target_link_libraries(DISTURBANCE DYNAMICS GLOBAL_ENVIRONMENT LOCAL_ENVIRONMENT ${S2E_LIBRARIES})
target_link_libraries(SIMULATION DYNAMICS GLOBAL_ENVIRONMENT LOCAL_ENVIRONMENT DISTURBANCE ${S2E_LIBRARIES})
target_link_libraries(GLOBAL_ENVIRONMENT ${CSPICE_LIB} ${S2E_LIBRARIES})
target_link_libraries(LOCAL_ENVIRONMENT GLOBAL_ENVIRONMENT ${CSPICE_LIB} ${S2E_LIBRARIES})
target_link_libraries(WRAPPER_NRLMSISE00 ${NRLMSISE00_LIB})

target_link_libraries(${PROJECT_NAME} DYNAMICS)
target_link_libraries(${PROJECT_NAME} DISTURBANCE)
target_link_libraries(${PROJECT_NAME} SIMULATION)
target_link_libraries(${PROJECT_NAME} GLOBAL_ENVIRONMENT LOCAL_ENVIRONMENT)
target_link_libraries(${PROJECT_NAME} RELATIVE_INFO)
target_link_libraries(${PROJECT_NAME} INI_ACC LOG_OUT SC_IO)
target_link_libraries(${PROJECT_NAME} COMPONENT)
target_link_libraries(${PROJECT_NAME} HILS_IO)

## C2A integration
if(USE_C2A)
  target_link_libraries(${PROJECT_NAME} C2A)
endif()

## HILS
if(USE_HILS)
  target_link_libraries(${PROJECT_NAME} ${WS2_32_LIB})
  set_target_properties(${PROJECT_NAME} PROPERTIES COMMON_LANGUAGE_RUNTIME "")
  set_target_properties(COMPONENT PROPERTIES COMMON_LANGUAGE_RUNTIME "")
  set_target_properties(DYNAMICS PROPERTIES COMMON_LANGUAGE_RUNTIME "")
  set_target_properties(DISTURBANCE PROPERTIES COMMON_LANGUAGE_RUNTIME "")
  set_target_properties(SIMULATION PROPERTIES COMMON_LANGUAGE_RUNTIME "")
  set_target_properties(GLOBAL_ENVIRONMENT PROPERTIES COMMON_LANGUAGE_RUNTIME "")
  set_target_properties(LOCAL_ENVIRONMENT PROPERTIES COMMON_LANGUAGE_RUNTIME "")
  set_target_properties(HILS_IO PROPERTIES COMMON_LANGUAGE_RUNTIME "")
  set_target_properties(RELATIVE_INFO PROPERTIES COMMON_LANGUAGE_RUNTIME "")
endif()

## Cmake debug
message("Cspice_LIB:  "${CSPICE_LIB})
message("nrlmsise00_LIB:  "${NRLMSISE00_LIB})

include(${S2E_CORE_DIR}/common.cmake)