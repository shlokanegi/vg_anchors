cmake_minimum_required(VERSION 3.10)
project(snarl_tree_builder)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -pthread")

# Path to libbdsg
set(LIBBDSG_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../libbdsg")

# Include directories
include_directories(
    ${LIBBDSG_DIR}/bdsg/include
    ${LIBBDSG_DIR}/bdsg/deps/libhandlegraph/src/include
    ${LIBBDSG_DIR}/bdsg/deps/sdsl-lite/include
    ${LIBBDSG_DIR}/bdsg/deps/sparsepp
)

# Link directories
link_directories(
    ${LIBBDSG_DIR}/lib
    ${LIBBDSG_DIR}
)

# Executable
add_executable(construct_snarl_tree_map construct_snarl_tree_map.cpp)

# Link libraries
target_link_libraries(construct_snarl_tree_map
    bdsg
    handlegraph
    sdsl
    divsufsort
    divsufsort64
    pthread
)

# Set RPATH so the executable can find shared libraries
set_target_properties(construct_snarl_tree_map PROPERTIES
    BUILD_RPATH "${LIBBDSG_DIR}:${LIBBDSG_DIR}/lib"
    INSTALL_RPATH "${LIBBDSG_DIR}:${LIBBDSG_DIR}/lib"
)

