# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.26

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/marius/PhD/CellMotility/agent_simulation

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/marius/PhD/CellMotility/agent_simulation/build

# Include any dependencies generated for this target.
include src/CMakeFiles/agent_simulation.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/CMakeFiles/agent_simulation.dir/compiler_depend.make

# Include the progress variables for this target.
include src/CMakeFiles/agent_simulation.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/agent_simulation.dir/flags.make

src/CMakeFiles/agent_simulation.dir/main.c.o: src/CMakeFiles/agent_simulation.dir/flags.make
src/CMakeFiles/agent_simulation.dir/main.c.o: /home/marius/PhD/CellMotility/agent_simulation/src/main.c
src/CMakeFiles/agent_simulation.dir/main.c.o: src/CMakeFiles/agent_simulation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/marius/PhD/CellMotility/agent_simulation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object src/CMakeFiles/agent_simulation.dir/main.c.o"
	cd /home/marius/PhD/CellMotility/agent_simulation/build/src && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT src/CMakeFiles/agent_simulation.dir/main.c.o -MF CMakeFiles/agent_simulation.dir/main.c.o.d -o CMakeFiles/agent_simulation.dir/main.c.o -c /home/marius/PhD/CellMotility/agent_simulation/src/main.c

src/CMakeFiles/agent_simulation.dir/main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/agent_simulation.dir/main.c.i"
	cd /home/marius/PhD/CellMotility/agent_simulation/build/src && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/marius/PhD/CellMotility/agent_simulation/src/main.c > CMakeFiles/agent_simulation.dir/main.c.i

src/CMakeFiles/agent_simulation.dir/main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/agent_simulation.dir/main.c.s"
	cd /home/marius/PhD/CellMotility/agent_simulation/build/src && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/marius/PhD/CellMotility/agent_simulation/src/main.c -o CMakeFiles/agent_simulation.dir/main.c.s

# Object files for target agent_simulation
agent_simulation_OBJECTS = \
"CMakeFiles/agent_simulation.dir/main.c.o"

# External object files for target agent_simulation
agent_simulation_EXTERNAL_OBJECTS =

src/agent_simulation: src/CMakeFiles/agent_simulation.dir/main.c.o
src/agent_simulation: src/CMakeFiles/agent_simulation.dir/build.make
src/agent_simulation: src/libcell_simulation_library.a
src/agent_simulation: src/libxoshiro_rng.a
src/agent_simulation: src/CMakeFiles/agent_simulation.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/marius/PhD/CellMotility/agent_simulation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable agent_simulation"
	cd /home/marius/PhD/CellMotility/agent_simulation/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/agent_simulation.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/agent_simulation.dir/build: src/agent_simulation
.PHONY : src/CMakeFiles/agent_simulation.dir/build

src/CMakeFiles/agent_simulation.dir/clean:
	cd /home/marius/PhD/CellMotility/agent_simulation/build/src && $(CMAKE_COMMAND) -P CMakeFiles/agent_simulation.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/agent_simulation.dir/clean

src/CMakeFiles/agent_simulation.dir/depend:
	cd /home/marius/PhD/CellMotility/agent_simulation/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/marius/PhD/CellMotility/agent_simulation /home/marius/PhD/CellMotility/agent_simulation/src /home/marius/PhD/CellMotility/agent_simulation/build /home/marius/PhD/CellMotility/agent_simulation/build/src /home/marius/PhD/CellMotility/agent_simulation/build/src/CMakeFiles/agent_simulation.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/agent_simulation.dir/depend

