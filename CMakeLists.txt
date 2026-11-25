cmake_minimum_required(VERSION 3.18)
project(helmer LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Dependencies: FLINT (>= 3.0), GMP, msolve's libneogb (>= 0.6).
find_library(FLINT_LIB flint REQUIRED)
find_library(GMP_LIB gmp REQUIRED)
find_library(NEOGB_LIB neogb REQUIRED)
find_path(FLINT_INCLUDE_DIR flint/flint.h REQUIRED)   # outside the compiler's default paths on Homebrew
find_path(GMP_INCLUDE_DIR gmp.h REQUIRED)
find_package(OpenMP)

add_library(helmer_core STATIC
  src/poly.cpp src/gb.cpp src/ideal.cpp src/primes.cpp src/whitney.cpp src/certify.cpp)
target_include_directories(helmer_core PUBLIC include ${FLINT_INCLUDE_DIR} ${GMP_INCLUDE_DIR})
target_link_libraries(helmer_core PUBLIC ${NEOGB_LIB} ${FLINT_LIB} ${GMP_LIB} m)
if(OpenMP_C_FOUND)
  target_link_libraries(helmer_core PUBLIC OpenMP::OpenMP_C)
endif()

option(HELMER_PYTHON "Build the Python extension" ON)
if(HELMER_PYTHON)
  find_package(Python 3.9 COMPONENTS Interpreter Development.Module REQUIRED)
  execute_process(COMMAND ${Python_EXECUTABLE} -m pybind11 --cmakedir
                  OUTPUT_VARIABLE pybind11_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
  find_package(pybind11 CONFIG REQUIRED)
  pybind11_add_module(_helmer src/python_bindings.cpp)
  target_link_libraries(_helmer PRIVATE helmer_core)
  install(TARGETS _helmer DESTINATION helmer)
endif()

option(HELMER_TESTS "Build C++ test drivers" ON)
if(HELMER_TESTS)
  foreach(t t_gb t_primes)
    add_executable(${t} tests/${t}.cpp)
    target_link_libraries(${t} PRIVATE helmer_core)
  endforeach()
endif()

option(HELMER_EXAMPLES "Build the C++ example (examples/cpp)" OFF)
if(HELMER_EXAMPLES)
  add_executable(stratify examples/cpp/stratify.cpp)
  target_link_libraries(stratify PRIVATE helmer_core)
endif()
