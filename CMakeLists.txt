cmake_minimum_required(VERSION 3.16)

# set the project name
project(KERO VERSION 0.1)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if (NOT CMAKE_CXX_FLAGS)
    set(CMAKE_CXX_FLAGS "-Wall -Ofast -g")
endif()

add_subdirectory("external/pthash")

add_library(kero
        src/kero_io.cpp
        src/util.cpp
        src/kero_mmap.cpp
)

add_custom_target(
        TurboPFor
        COMMAND make -C "${PROJECT_SOURCE_DIR}/external/TurboPFor"
        COMMENT "Building TurboPFor..."
        VERBATIM
)
add_dependencies(kero TurboPFor)

target_include_directories(kero
        PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/external/pthash/include>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/external/logfault>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/external/TurboPFor/include>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
)

set(TURBO_PFOR_LIB
        ${PROJECT_SOURCE_DIR}/external/TurboPFor/libic.a
)

target_link_libraries(kero
        PTHASH
        ${TURBO_PFOR_LIB}
)