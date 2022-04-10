add_test( HelloTest.SameStrings /home/marius/PhD/CellMotility/agent_simulation/build/test/run_unit_tests [==[--gtest_filter=HelloTest.SameStrings]==] --gtest_also_run_disabled_tests)
set_tests_properties( HelloTest.SameStrings PROPERTIES WORKING_DIRECTORY /home/marius/PhD/CellMotility/agent_simulation/build/test SKIP_REGULAR_EXPRESSION [==[\[  SKIPPED \]]==])
set( run_unit_tests_TESTS HelloTest.SameStrings)
