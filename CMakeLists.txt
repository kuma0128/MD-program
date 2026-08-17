cmake_minimum_required(VERSION 3.16)
project(md_program LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(md STATIC
    EWALD/ewald.cpp
    MC/monte_carlo.cpp
    RATTLE/rattle.cpp
    ROTATE/rotate.cpp
    SHAKE/shake.cpp
)
target_include_directories(md PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(md PRIVATE -Wall -Wextra -Wpedantic)

add_executable(gaussian GAUSS/main.cpp)
target_include_directories(gaussian PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(gaussian PRIVATE -Wall -Wextra -Wpedantic)

enable_testing()
add_executable(md_tests tests/test_md.cpp)
target_link_libraries(md_tests PRIVATE md)
target_compile_options(md_tests PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME md_tests COMMAND md_tests)
