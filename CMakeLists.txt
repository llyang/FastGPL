cmake_minimum_required(VERSION 3.16)
project(FastGPL VERSION 0.1)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

include_directories(include)

add_subdirectory(src)

add_subdirectory(test)


