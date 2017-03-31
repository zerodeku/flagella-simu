cmake_minimum_required(VERSION 3.5)
project(flagella_simu)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

message("${CMAKE_CXX_FLAGS}")

set(flagella_simu_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
include_directories(${flagella_simu_INCLUDE_DIR})
include_directories(armadillo/include)

file(GLOB flagella_SRC "${PROJECT_SOURCE_DIR}/src/*.cpp")

add_subdirectory(debug)
add_subdirectory(src)

set(${ARMADILLO_LIBRARIES} armadillo/lib/libarmadillo.dll.a
        armadillo/lib/libarmadillo.dll
        )

target_link_libraries(simulation ${ARMADILLO_LIBRARIES})