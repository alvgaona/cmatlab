cmake_minimum_required(VERSION 3.10)

project(some)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(OpenCV REQUIRED)
find_package(OpenMP REQUIRED)

include_directories(${OpenCV_INCLUDE_DIRS})

# Add executable target with source files listed in the SOURCES variable
file(GLOB SOURCES "src/resizeim/*.cpp")
add_executable(${PROJECT_NAME} ${SOURCES} "src/main.cpp")

target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)

# Link all .so files in the lib directory
file(GLOB SO_LIBS "lib/*.so*")
target_link_libraries(${PROJECT_NAME} ${OpenCV_LIBS} ${SO_LIBS})
