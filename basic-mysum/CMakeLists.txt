cmake_minimum_required(VERSION 3.10)

project(mysum)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(opencv4 REQUIRED)

include_directories(include)

# Add executable target with source files listed in the SOURCES variable
file(GLOB SOURCES "src/*.cpp")
add_executable(${PROJECT_NAME} ${SOURCES})
