# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.23

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
include test/CMakeFiles/run_unit_tests.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include test/CMakeFiles/run_unit_tests.dir/compiler_depend.make

# Include the progress variables for this target.
include test/CMakeFiles/run_unit_tests.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/run_unit_tests.dir/flags.make

test/CMakeFiles/run_unit_tests.dir/unit_tests.cc.o: test/CMakeFiles/run_unit_tests.dir/flags.make
test/CMakeFiles/run_unit_tests.dir/unit_tests.cc.o: ../test/unit_tests.cc
test/CMakeFiles/run_unit_tests.dir/unit_tests.cc.o: test/CMakeFiles/run_unit_tests.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/marius/PhD/CellMotility/agent_simulation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/run_unit_tests.dir/unit_tests.cc.o"
	cd /home/marius/PhD/CellMotility/agent_simulation/build/test && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/CMakeFiles/run_unit_tests.dir/unit_tests.cc.o -MF CMakeFiles/run_unit_tests.dir/unit_tests.cc.o.d -o CMakeFiles/run_unit_tests.dir/unit_tests.cc.o -c /home/marius/PhD/CellMotility/agent_simulation/test/unit_tests.cc

test/CMakeFiles/run_unit_tests.dir/unit_tests.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/run_unit_tests.dir/unit_tests.cc.i"
	cd /home/marius/PhD/CellMotility/agent_simulation/build/test && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/marius/PhD/CellMotility/agent_simulation/test/unit_tests.cc > CMakeFiles/run_unit_tests.dir/unit_tests.cc.i

test/CMakeFiles/run_unit_tests.dir/unit_tests.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/run_unit_tests.dir/unit_tests.cc.s"
	cd /home/marius/PhD/CellMotility/agent_simulation/build/test && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/marius/PhD/CellMotility/agent_simulation/test/unit_tests.cc -o CMakeFiles/run_unit_tests.dir/unit_tests.cc.s

# Object files for target run_unit_tests
run_unit_tests_OBJECTS = \
"CMakeFiles/run_unit_tests.dir/unit_tests.cc.o"

# External object files for target run_unit_tests
run_unit_tests_EXTERNAL_OBJECTS =

test/run_unit_tests: test/CMakeFiles/run_unit_tests.dir/unit_tests.cc.o
test/run_unit_tests: test/CMakeFiles/run_unit_tests.dir/build.make
test/run_unit_tests: src/libcell_simulation_library.a
test/run_unit_tests: test/CMakeFiles/run_unit_tests.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/marius/PhD/CellMotility/agent_simulation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable run_unit_tests"
	cd /home/marius/PhD/CellMotility/agent_simulation/build/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/run_unit_tests.dir/link.txt --verbose=$(VERBOSE)
	cd /home/marius/PhD/CellMotility/agent_simulation/build/test && /usr/bin/cmake -D TEST_TARGET=run_unit_tests -D TEST_EXECUTABLE=/home/marius/PhD/CellMotility/agent_simulation/build/test/run_unit_tests -D TEST_EXECUTOR= -D TEST_WORKING_DIR=/home/marius/PhD/CellMotility/agent_simulation/build/test -D TEST_EXTRA_ARGS= -D TEST_PROPERTIES= -D TEST_PREFIX= -D TEST_SUFFIX= -D TEST_FILTER= -D NO_PRETTY_TYPES=FALSE -D NO_PRETTY_VALUES=FALSE -D TEST_LIST=run_unit_tests_TESTS -D CTEST_FILE=/home/marius/PhD/CellMotility/agent_simulation/build/test/run_unit_tests[1]_tests.cmake -D TEST_DISCOVERY_TIMEOUT=5 -D TEST_XML_OUTPUT_DIR= -P /usr/share/cmake/Modules/GoogleTestAddTests.cmake

# Rule to build all files generated by this target.
test/CMakeFiles/run_unit_tests.dir/build: test/run_unit_tests
.PHONY : test/CMakeFiles/run_unit_tests.dir/build

test/CMakeFiles/run_unit_tests.dir/clean:
	cd /home/marius/PhD/CellMotility/agent_simulation/build/test && $(CMAKE_COMMAND) -P CMakeFiles/run_unit_tests.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/run_unit_tests.dir/clean

test/CMakeFiles/run_unit_tests.dir/depend:
	cd /home/marius/PhD/CellMotility/agent_simulation/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/marius/PhD/CellMotility/agent_simulation /home/marius/PhD/CellMotility/agent_simulation/test /home/marius/PhD/CellMotility/agent_simulation/build /home/marius/PhD/CellMotility/agent_simulation/build/test /home/marius/PhD/CellMotility/agent_simulation/build/test/CMakeFiles/run_unit_tests.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/run_unit_tests.dir/depend

