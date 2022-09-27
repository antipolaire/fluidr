#ifndef ALEPH_TESTS_BASE_HH__
#define ALEPH_TESTS_BASE_HH__

#define CMAKE_SOURCE_DIR         "/home/david/CLionProjects/fluid_solver_raw"
#define CMAKE_CURRENT_BINARY_DIR "/home/david/CLionProjects/fluid_solver_raw/cmake-build-debug/test"
/* #undef CMAKE_PROJECT_BINARY_DIR */

#define CMAKE_MAJOR_VERSION      "3"
#define CMAKE_MINOR_VERSION      "23"
#define CMAKE_PATCH_VERSION      "2"

#include <iostream>
#include <stdexcept>
#include <string>

namespace aleph
{

#define ALEPH_ASSERT_EQUAL( x, y )                                  \
{                                                                   \
  if( ( x ) != ( y ) )                                              \
  {                                                                 \
    throw std::runtime_error(   std::string( __FILE__ )             \
                              + std::string( ":" )                  \
                              + std::to_string( __LINE__ )          \
                              + std::string( " in " )               \
                              + std::string( __PRETTY_FUNCTION__ )  \
                              + std::string( ": " )                 \
                              + std::to_string( ( x ) )             \
                              + std::string( " != " )               \
                              + std::to_string( ( y ) )             \
    );                                                              \
  }                                                                 \
}


#define ALEPH_ASSERT_THROW( condition )                             \
{                                                                   \
  if( !( condition ) )                                              \
  {                                                                 \
    throw std::runtime_error(   std::string( __FILE__ )             \
                              + std::string( ":" )                  \
                              + std::to_string( __LINE__ )          \
                              + std::string( " in " )               \
                              + std::string( __PRETTY_FUNCTION__ )  \
    );                                                              \
  }                                                                 \
}

#define ALEPH_EXPECT_EXCEPTION( expression, exception )             \
{                                                                   \
  try                                                               \
  {                                                                 \
    ( expression );                                                 \
  }                                                                 \
  catch( exception& e )                                             \
  {                                                                 \
  }                                                                 \
  catch( ... )                                                      \
  {                                                                 \
    throw std::runtime_error(   std::string( __FILE__ )             \
                              + std::string( ":" )                  \
                              + std::to_string( __LINE__ )          \
                              + std::string( " in " )               \
                              + std::string( __PRETTY_FUNCTION__ )  \
    );                                                              \
  }                                                                 \
}

#define ALEPH_TEST_BEGIN( name )\
{\
  std::cerr << "-- Running test \"" << name << "\"...";\
}

#define ALEPH_TEST_END() \
{\
  std::cerr << "finished\n";\
}


}

#endif
