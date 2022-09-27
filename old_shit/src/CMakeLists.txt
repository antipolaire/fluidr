
# ####
# Package: Eigen3
# ####
find_package(Eigen3 3.3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})
link_libraries(Eigen3::Eigen)

# ####
# Package: SDL2
# ####
find_package(SDL2 REQUIRED)
include_directories(${SDL2_INCLUDE_DIRS})
link_libraries(${SDL2_LIBRARY} ${SDLIMAGE_LIBRARY})

# ####
# Package: NetPBM
# ####
find_package(PBM REQUIRED)
include_directories(${NETPBM_INCLUDE_DIR})
link_libraries(${NETPBM_LIBRARY})

# ####
# Cinder OpenGL Wrapper
# ####

# find_package(CINDER REQUIRED)
# include_directories(${CINDER_INCLUDE_DIRS})
# link_libraries(${CINDER_LIBRARY})

# ####
# Suite Parse containing the awesome Cholmod solver
# ####

# include_directories(/Users/david.baldin/projects/tssrio/SuiteSparse/CHOLMOD/Include)
# link_libraries(/Users/david.baldin/projects/tssrio/SuiteSparse/CHOLMOD/Lib/libcholmod.a)
# include_directories("/usr/include/suitesparse")
# link_directories(/usr/local/lib)
#link_libraries(libcholmod.a)


## PATH SETTINGS

get_filename_component(APP_PATH "${PROJECT_SOURCE_DIR}" ABSOLUTE)
get_filename_component(CINDER_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../Cinder" ABSOLUTE)

## INCLUDE CINDER
include("${CINDER_PATH}/proj/cmake/modules/cinderMakeApp.cmake")

ci_make_app(
        APP_NAME "solver_cinder"
        CINDER_PATH ${CINDER_PATH}
        SOURCES
        "${APP_PATH}/src/main_cinder.cpp"
        "${APP_PATH}/src/field2d_v2.hpp"
        "${APP_PATH}/src/solver.hpp"
        "${APP_PATH}/src/solver_v2.hpp"
        INCLUDES "${APP_PATH}/include"
        RESOURCES
        "${APP_PATH}/resource/assets/shaders/common/passthr.vert"
        "${APP_PATH}/resource/assets/shaders/common/passthr.frag"
        "${APP_PATH}/resource/assets/shaders/common/dots.geom"
        LIBRARIES Eigen3::Eigen)