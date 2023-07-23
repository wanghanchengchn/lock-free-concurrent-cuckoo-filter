cmake_minimum_required(VERSION 3.5)
project(concurrent-filter LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
# Write compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Threads)

# add googletest
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.13.0.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -mavx2 -mbmi -mbmi2") # avx2, bmi, bmi2 for vqf

# string (TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UC)
# if (CMAKE_BUILD_TYPE_UC STREQUAL "RELEASE")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
# endif()

add_subdirectory(src)
add_subdirectory(test)