#Unit testing
enable_testing()

##CMake native testing
#add_test(NAME Runs COMMAND agent_simulation "/home/marius/PhD/CellMotility/agent_simulation/output/test/test_parameters")

#set_tests_properties(Runs PROPERTIES PASS_REGULAR_EXPRESSION "0")

#Google test testing
#Download gtest (should probably link to newer version) URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
include(FetchContent)
FetchContent_Declare(
    googletest URL https://github.com/google/googletest/archive/master.zip
    
)


add_executable(run_unit_tests unit_tests.cc)


#gtest_main makes sure you don't need your own main function to run tests
target_link_libraries(run_unit_tests PUBLIC
    gtest_main
    gtest 
    cell_simulation_library
    )

include(GoogleTest)
gtest_discover_tests(run_unit_tests) 
