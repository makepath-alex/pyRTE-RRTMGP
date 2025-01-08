cmake_minimum_required(VERSION 3.18...3.26)
project(rte-rrtmgp)

add_custom_target(compile-rte-rrtmgp ALL)

# Determine platform
set(LINUX FALSE)
if (UNIX AND NOT APPLE)
    set(LINUX TRUE)
endif()

# Compile Python bindings with pybind11
find_package(pybind11 REQUIRED)

set(TARGET_NAME pyrte_rrtmgp)
set(OUTPUT_DIRECTORY pyrte_rrtmgp)
set(SOURCES ${CMAKE_SOURCE_DIR}/pybind_interface.cpp)

pybind11_add_module(${TARGET_NAME} ${SOURCES})

target_include_directories(${TARGET_NAME} PUBLIC
    $ENV{CONDA_PREFIX}/include
)

target_link_directories(${TARGET_NAME} PUBLIC
    $ENV{CONDA_PREFIX}/lib
)

target_link_libraries(${TARGET_NAME} PUBLIC
    rte
    rrtmgp
)

target_compile_definitions(${TARGET_NAME} PRIVATE
    VERSION_INFO=${VERSION_INFO}
    DBL_EPSILON=5.8e-2
    DCMAKE_LIBRARY_OUTPUT_DIRECTORY=${OUTPUT_DIRECTORY}
)

# Add these checks after the initial Linux check
set(APPLE_ARM FALSE)
if(APPLE)
    if(CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64")
        set(APPLE_ARM TRUE)
    endif()
endif()

if (${LINUX})
    target_link_libraries(${TARGET_NAME} PUBLIC gfortran)
elseif(APPLE)
    # On macOS, explicitly link against gfortran runtime
    if(APPLE_ARM)
        target_link_directories(${TARGET_NAME} PUBLIC /opt/homebrew/lib/gcc/current)
    else()
        target_link_directories(${TARGET_NAME} PUBLIC /usr/local/lib/gcc/current)
    endif()
    target_link_libraries(${TARGET_NAME} PUBLIC gfortran)
endif()

# The install directory is the output (wheel) directory
install(TARGETS ${TARGET_NAME} DESTINATION ${OUTPUT_DIRECTORY})
