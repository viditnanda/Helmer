# core_gba exposes internal structs, so compile the adapter and backend
# against the same source. The default (non-trace) build uses system neogb.
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)
endif()
include(FetchContent)
FetchContent_Declare(helmer_msolve
  URL https://codeload.github.com/algebraic-solving/msolve/tar.gz/0ae16f24d3c6ef2457448e39f2de33a5ac1d6b72
  URL_HASH SHA256=8403ca3cfa0063c88e02872b9efffdf1053cb3f741a8abd454a0ef58a34065ae)
FetchContent_MakeAvailable(helmer_msolve)
set(neogb_src "${helmer_msolve_SOURCE_DIR}/src/neogb")
set(neogb_build "${CMAKE_CURRENT_BINARY_DIR}/trace-neogb")

# Patch a private copy of the amalgamation's affected includes. Never
# modify a caller's offline source tree. Restore pristine inputs on reconfigure.
foreach(source gb.c f4.c hash.c)
  configure_file("${neogb_src}/${source}" "${neogb_build}/${source}" COPYONLY)
endforeach()
find_package(Git REQUIRED)
set(trace_patch "${CMAKE_CURRENT_LIST_DIR}/msolve-trace-cleanup.patch")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${trace_patch}")
execute_process(
  COMMAND "${CMAKE_COMMAND}" -E env "GIT_CEILING_DIRECTORIES=${CMAKE_CURRENT_BINARY_DIR}"
          "${GIT_EXECUTABLE}" apply --no-index -p3 "${trace_patch}"
  WORKING_DIRECTORY "${neogb_build}" RESULT_VARIABLE patch_status)
if(NOT patch_status EQUAL 0)
  message(FATAL_ERROR "Could not apply the pinned msolve replay fixes")
endif()

add_library(helmer_trace_backend STATIC "${neogb_build}/gb.c")
set_target_properties(helmer_trace_backend PROPERTIES C_VISIBILITY_PRESET hidden C_STANDARD 11 C_STANDARD_REQUIRED YES)
target_include_directories(helmer_trace_backend PUBLIC "${neogb_src}")
target_link_libraries(helmer_trace_backend PRIVATE ${GMP_LIB} m)
if(OpenMP_C_FOUND)
  target_link_libraries(helmer_trace_backend PRIVATE OpenMP::OpenMP_C)
endif()
