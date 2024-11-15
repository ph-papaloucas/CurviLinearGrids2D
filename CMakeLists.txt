cmake_minimum_required(VERSION 3.10)
project(CurviLinearGrids2D)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Only main.cpp is the executable source
add_executable(${PROJECT_NAME} main.cpp)

# Include directories for header files
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Link the necessary source files (Grid2D.cpp) to the executable
target_sources(${PROJECT_NAME} PRIVATE Grid2D.cpp utils.cpp)


add_executable(grid2 grid2.cpp)