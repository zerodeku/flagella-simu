cmake_minimum_required(VERSION 3.5)
set(SOURCE_FILES sim_driver.cpp)
add_executable(simulation ${flagella_SRC} ${SOURCE_FILES})

find_package(Boost COMPONENTS thread date_time program_options filesystem system REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

if (WIN32 AND __COMPILER_GNU)
    # mingw-gcc fails to link boost::thread
    add_definitions(-DBOOST_THREAD_USE_LIB)
endif (WIN32 AND __COMPILER_GNU)

target_link_libraries(simulation
        ${CMAKE_THREAD_LIBS_INIT}
        ${Boost_LIBRARIES}
        )