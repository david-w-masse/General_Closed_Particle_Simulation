﻿# CMakeList.txt : CMake project for General_Closed_Particle_Simulation, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("General_Closed_Particle_Simulation")

# Add source to this project's executable.
add_executable (General_Closed_Particle_Simulation "General_Closed_Particle_Simulation.cpp" "General_Closed_Particle_Simulation.h" "Initialization_2D.cpp" "Data_Loader.cpp" "Initialization_3D.cpp" "Calculate_Forces_2D.cpp" "Calculate_Forces_3D.cpp")

# let the compiler use avx instructions
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET General_Closed_Particle_Simulation PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
