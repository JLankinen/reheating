cmake_minimum_required(VERSION 3.16)
project(KinationReheatingScenario)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(include)

file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS src/*.cpp)

add_executable(reheating ${SOURCES})

find_package(Boost REQUIRED)

if(Boost_FOUND)
    message(STATUS "Boost found at ${Boost_INCLUDE_DIRS}")
    target_include_directories(reheating PRIVATE ${Boost_INCLUDE_DIRS})
    target_link_libraries(reheating PRIVATE ${Boost_LIBRARIES})
endif()

# Warnings
target_compile_options(reheating PRIVATE -Wall -Wextra -Wpedantic)