cmake_minimum_required(VERSION 3.10)
project(fftSqueeze)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_COMPILER clang++)

include_directories(include)

file(GLOB SOURCES "src/*.cpp")

add_executable(myprogram ${SOURCES})
