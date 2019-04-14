cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(ColorConverter)

#option(WITH_TESTS "Whether tests should be built and run" ON)

# build the main target
add_executable(${PROJECT_NAME}
  example/main.cpp
)

target_include_directories(${PROJECT_NAME}
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/include
)

set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "-Wall -Weverything -Wno-c++98-compat -Wno-padded -Wno-documentation")
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD_REQUIRED ON)
