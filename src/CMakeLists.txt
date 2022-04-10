# Create libraries
add_library(cell_simulation_library cell_simulation.c)
add_library(xoshiro_rng xoshiro_rng.c)

#add the executable
add_executable(agent_simulation main.c)

# add compile options: math, gsl, gslcblas (for optimisation routines)
target_link_libraries(agent_simulation PUBLIC 
                        cell_simulation_library 
                        xoshiro_rng
                        m 
                        gsl 
                        gslcblas
                        )

# add the binary tree to the search path for the include files so 
# that we will find agent_simulation_config.h
#target_include_directories("${PROJECT_SOURCE_DIR}/src")

#target_include_directories(agent_simulation PUBLIC
                            #"${PROJECT_BINARY_DIR}"
                            #)
