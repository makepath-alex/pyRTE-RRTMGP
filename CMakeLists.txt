cmake_minimum_required(VERSION 3.18...3.26)
project(rte-rrtmgp)

# Helper function to check if a value is in a list
function(check_valid_option variable_name value valid_options)
    list(FIND ${valid_options} ${value} index)
    if (index EQUAL -1)
        message(FATAL_ERROR
            "${variable_name} has an invalid value: '${value}'.\n"
            "Valid options are: ${${valid_options}}.")
    endif()
endfunction()

set(VALID_RTE_KERNELS "ACCEL" "EXTERN" "DEFAULT")
set(RTE_KERNELS "DEFAULT")

# Set the RTE_KERNELS value to upper case only
string(TOUPPER "${RTE_KERNELS}" RTE_KERNELS)

# Validate RTE_KERNELS
check_valid_option("RTE_KERNELS" "${RTE_KERNELS}" "VALID_RTE_KERNELS")

message(STATUS "RTE_KERNELS: ${RTE_KERNELS}")

# Python paths
set(PYTHON_EXECUTABLE "/home/vkm/code/makepath/pyRTE-RRTMGP/.venv/bin/python")
set(PYTHON_INCLUDE_DIR "/home/vkm/code/makepath/pyRTE-RRTMGP/.venv/include/python3.10")
set(PYTHON_LIBRARY "/home/vkm/code/makepath/pyRTE-RRTMGP/.venv/lib/libpython3.10.so")

set(RTE_RRTMGP_DIR  ${CMAKE_SOURCE_DIR}/build/rte-rrtmgp/)
set(RRTMGP_LIB_PATH ${RTE_RRTMGP_DIR}/build/rrtmgp-frontend/librrtmgp.a)
set(RTE_LIB_PATH    ${RTE_RRTMGP_DIR}/build/rte-frontend/librte.a)

# Ensure the libraries are built before linking
add_custom_command(
    OUTPUT ${RRTMGP_LIB_PATH} ${RTE_LIB_PATH}
    COMMAND ${CMAKE_COMMAND} -E echo "Waiting for rte-rrtmgp libraries to be built..."
    DEPENDS rte-rrtmgp
)

add_custom_target(wait_for_rte_rrtmgp_libs ALL
    DEPENDS ${RRTMGP_LIB_PATH} ${RTE_LIB_PATH}
)

# Custom target for always running
add_custom_target(compile-rte-rrtmgp ALL)

# Determine platform
set(LINUX FALSE)
if (UNIX AND NOT APPLE)
    set(LINUX TRUE)
endif()

# Fortran compiler selection
set(PREFERRED_FC_COMPILERS gfortran ifort gfortran-10 gfortran-11 gfortran-12 f77)

foreach(compiler IN LISTS PREFERRED_FC_COMPILERS)
    find_program(FOUND_COMPILER NAMES ${compiler})
    if(FOUND_COMPILER)
        set(CMAKE_Fortran_COMPILER ${FOUND_COMPILER})
        message(STATUS "Using Fortran compiler: ${FOUND_COMPILER}")
        break()
    endif()
endforeach()

if(NOT CMAKE_Fortran_COMPILER)
    message(FATAL_ERROR "No suitable Fortran compiler found")
endif()

# Include ExternalProject module
include(ExternalProject REQUIRED)

# Build command for external project
set(BUILD_COMMAND_STRING
    "rm -rf build && "
    "mkdir build && "
    "cd build && cmake .. -DRTE_KERNELS=${RTE_KERNELS} && "
    "cmake --build . -j8"
)

ExternalProject_Add(
    rte-rrtmgp
    GIT_REPOSITORY /home/vkm/code/makepath/rte-rrtmgp/.git
    GIT_TAG origin/cmake
    GIT_SHALLOW TRUE
    SOURCE_DIR ${RTE_RRTMGP_DIR}
    BUILD_IN_SOURCE TRUE
    CONFIGURE_COMMAND ""
    BUILD_COMMAND eval ${BUILD_COMMAND_STRING}
    INSTALL_COMMAND ""
)

# Compile Python bindings with pybind11
find_package(pybind11 REQUIRED)

set(TARGET_NAME pyrte_rrtmgp)
set(SOURCES ${CMAKE_SOURCE_DIR}/pybind_interface.cpp)

pybind11_add_module(${TARGET_NAME} ${SOURCES})

add_dependencies(${TARGET_NAME} wait_for_rte_rrtmgp_libs)

target_include_directories(${TARGET_NAME} PUBLIC
    ${RTE_RRTMGP_DIR}/rte-kernels/api/
    ${RTE_RRTMGP_DIR}/rrtmgp-kernels/api/
    ${CMAKE_Fortran_MODULE_DIRECTORY}
    ${RTE_RRTMGP_DIR}/build/modules
)

target_link_directories(${TARGET_NAME} PUBLIC
    ${RTE_RRTMGP_DIR}/build
    ${RTE_RRTMGP_DIR}/build/rte-frontend
    ${RTE_RRTMGP_DIR}/build/rrtmgp-frontend
    ${RTE_RRTMGP_DIR}/build/rte-kernels
    ${RTE_RRTMGP_DIR}/build/rrtmgp-kernels
)

target_link_libraries(${TARGET_NAME} PUBLIC
    librtekernels.a
    librte.a
    librrtmgpkernels.a
    librrtmgp.a
)

target_compile_definitions(${TARGET_NAME} PRIVATE
    VERSION_INFO=${VERSION_INFO}
    DBL_EPSILON=5.8e-2
    DCMAKE_LIBRARY_OUTPUT_DIRECTORY=pyrte_rrtmgp
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
install(TARGETS pyrte_rrtmgp DESTINATION pyrte_rrtmgp)
