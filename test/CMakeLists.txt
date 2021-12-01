
add_executable(arc_test arc_test.c)
add_executable(interactive_arc_test interactive_arc_test.c)
add_executable(interactive_arc_front_test interactive_arc_front_test.c)
add_executable(openmp_test openmp_test.c)
target_link_libraries(arc_test PRIVATE arc)
target_link_libraries(interactive_arc_test PRIVATE arc)
target_link_libraries(interactive_arc_front_test  PRIVATE arc)
target_link_libraries(openmp_test PRIVATE arc)

add_test(arc_test arc_test)
add_test(openmp_test openmp_test)
