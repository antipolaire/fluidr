CONFIGURE_FILE(Base.hh.in ${CMAKE_SOURCE_DIR}/test/Base.hh)

include_directories(${TEST_SOURCE_DIR}/src)

# ####
# Eigen library for vector algebra and alike
# ####

find_package(Eigen3 3.3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})
link_libraries(Eigen3::Eigen)

# ####
# Suite Parse containing the awesome Cholmod solver
# ####

# include_directories("/usr/local/Cellar/suite-sparse/5.3.0_1/include")
# link_directories("/usr/local/Cellar/suite-sparse/5.3.0_1/lib")
include_directories("/usr/include/suitesparse")
# link_directories(/usr/local/lib)

# ####
# Add the main test file to the executable
# ####

add_executable(test_dummy test_dummy.cpp)
add_executable(test_field2d test_field2d.cpp)
add_executable(test_field2d_v2 test_field2d_v2.cpp)
add_executable(test_all_solver_performance test_all_solver_performance.cpp)
add_executable(test_mac_grid test_mac_grid.cpp)
add_executable(test_vector_translation test_vector_translation.cpp)
add_executable(test_flip_reference test_flip_reference.cpp)
add_executable(test_interpolation test_interpolation.cpp)
add_executable(test_integration test_integration.cpp)
add_executable(test_modular_solver test_modular_solver.cpp)
add_executable(test_gfx test_gfx.cpp)
add_executable(test_solver test_solver.cpp)
add_executable(test_solver_v2 test_solver_v2.cpp)
add_executable(test_eigen_matrix test_eigen_matrix.cpp)

# ####
# Link the cholmod library from suite parse to our executable
# ####
# link_libraries(libcholmod.a)

# target_link_libraries(test_solver libcholmod.a)

add_test(NAME dummy COMMAND test_dummy)
add_test(NAME field2d COMMAND test_field2d)
add_test(NAME test_field2d_v2 COMMAND test_field2d_v2)
add_test(NAME all_solver_performance COMMAND test_all_solver_performance)
add_test(NAME mac_grid COMMAND test_mac_grid)
add_test(NAME vector_translation COMMAND test_vector_translation)
add_test(NAME test_flip_reference COMMAND test_flip_reference)
add_test(NAME test_interpolation COMMAND test_interpolation)
add_test(NAME test_integration COMMAND test_integration)
add_test(NAME test_modular_solver COMMAND test_modular_solver)
add_test(NAME test_gfx COMMAND test_gfx)
add_test(NAME test_solver COMMAND test_solver)
add_test(NAME test_solver_v2 COMMAND test_solver_v2)
add_test(NAME test_eigen_matrix COMMAND test_eigen_matrix)