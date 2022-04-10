cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(agent_simulation VERSION 1.0)

# configure a header file to pass some of the CMake settings
# to the source code
#configure_file(agent_simulation_config.h.in agent_simulation_config.h)

#Set GCC version 11, since gtest needs that
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")

#add the main directory
add_subdirectory(src)

#add the test directory
add_subdirectory(test)



          

                            
