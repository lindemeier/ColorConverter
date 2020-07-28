cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(ColorConverter)

add_library(${PROJECT_NAME} INTERFACE)
target_sources(${PROJECT_NAME}
  INTERFACE ${PROJECT_SOURCE_DIR}/include/ColorConverter/ColorConverter.hxx
)

target_include_directories(${PROJECT_NAME}
  INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  INTERFACE $<INSTALL_INTERFACE:include/>
)


#option(WITH_TESTS "Whether tests should be built and run" ON)

# build the main target
add_executable(${PROJECT_NAME}_example
  example/main.cpp
)

target_link_libraries(${PROJECT_NAME}_example
  ${PROJECT_NAME}
)

add_dependencies(${PROJECT_NAME}_example
  ${PROJECT_NAME}
)


set_target_properties(${PROJECT_NAME}_example PROPERTIES COMPILE_FLAGS "-Wall -Weverything -Wno-c++98-compat -Wno-padded -Wno-documentation")
set_target_properties(${PROJECT_NAME}_example PROPERTIES CXX_STANDARD 17)
set_target_properties(${PROJECT_NAME}_example PROPERTIES CXX_STANDARD_REQUIRED ON)
