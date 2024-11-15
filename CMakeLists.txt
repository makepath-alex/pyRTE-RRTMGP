cmake_minimum_required(VERSION 3.18...3.26)
project(rte-rrtmgp)

SET(RTE_RRTMGP_PATH /home/vkm/code/makepath/rte-rrtmgp/.git) # https://github.com/earth-system-radiation/rte-rrtmgp.git
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

# Include ExternalProject module
include(ExternalProject REQUIRED)

# Build command for external project
set(BUILD_COMMAND_STRING
    "rm -rf build && "
    "mkdir build && "
    "cd build && cmake .. -DRTE_KERNELS=default && "
    "cmake --build . -j8"
    # "make -C build -j ${N}"
)

ExternalProject_Add(
    rte-rrtmgp
    GIT_REPOSITORY ${RTE_RRTMGP_PATH}
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
)

target_link_directories(${TARGET_NAME} PUBLIC
    ${RTE_RRTMGP_DIR}/build
    ${RTE_RRTMGP_DIR}/build/rte-frontend
    ${RTE_RRTMGP_DIR}/build/rrtmgp-frontend
    ${RTE_RRTMGP_DIR}/build/rte-kernels
    ${RTE_RRTMGP_DIR}/build/rrtmgp-kernels
)

target_link_libraries(${TARGET_NAME} PUBLIC
    rte
    rrtmgp
    rtekernels
    rrtmgpkernels
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
