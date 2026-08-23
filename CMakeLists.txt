cmake_minimum_required(VERSION 3.12)
project(my_simulation CXX C) # Put your project name and languages here (C++ == CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CODE_DIR ${CMAKE_SOURCE_DIR}/C)

include(ExternalProject) # Include the CMake command ExternalProject

# Search for system packages
find_package(Eigen3 3.4 REQUIRED) # Make sure that the compiler goes and finds Eigen with min version 3.4 in our system packages, define the path variable ${EIGEN3_INCLUDE_DIR}
find_package(lcm REQUIRED)
find_package(netCDF CONFIG REQUIRED)


# Add include directories, these are the directories when you #include "package"
include_directories(
    ${EIGEN3_INCLUDE_DIR}
    ${CODE_DIR}
    ${CODE_DIR}/include
    ${CODE_DIR}/include/transforms
    ${CODE_DIR}/lib
    ${CMAKE_SOURCE_DIR}/LCM
)

# Add the correct locations for GeographicLib to be found by the linker
# list(APPEND CMAKE_PREFIX_PATH
#     "${CODE_DIR}/lib/geographiclib"
# )
# find_package(GeographicLib REQUIRED)

# ----------------------------------------------------------------------------------------------------------
# Add the Spice kernels package, recompile if the CPU type has changed, otherwise do not recompile
execute_process(
    COMMAND uname -m
    OUTPUT_VARIABLE HOST_ARCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

if(NOT HOST_ARCH MATCHES "x86_64")
    message(STATUS "Non-x86 architecture (${HOST_ARCH}), patching CSPICE to remove -m64")
    file(WRITE "${CMAKE_BINARY_DIR}/patch_cspice.sh" "
#!/bin/sh
find '${CODE_DIR}/lib/spice/cspice' -name 'mkprodct.csh' | while read f; do
    sed -i 's/-m64//g' \"$f\"
done
")
    set(CSPICE_PATCH_COMMAND sh "${CMAKE_BINARY_DIR}/patch_cspice.sh")
else()
    message(STATUS "x86_64 detected, no CSPICE patch needed")
    set(CSPICE_PATCH_COMMAND ${CMAKE_COMMAND} -E echo "Skipping patch")  # proper no-op
endif()
# ----------------------------------------------------------------------------------------------------------

# Now that CSpice has been compiled it can be added as a library and compiled using its own make files
ExternalProject_Add(
    cspice_build
    SOURCE_DIR "${CODE_DIR}/lib/spice/cspice"
    CONFIGURE_COMMAND ""
    PATCH_COMMAND ${CSPICE_PATCH_COMMAND}
    BUILD_COMMAND ${CSH_EXECUTABLE} "${CODE_DIR}/lib/spice/cspice/makeall.csh"
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND ""
)
set(CSPICE_LIB "${CODE_DIR}/lib/spice/cspice/lib/cspice.a")


# We can also include an external project and use its own compile scripts with FetchContent
include(FetchContent)

# GeographicLib reads these; set before it's pulled in. Same two options
# the run script already passed, so we know they're the right names.
set(BUILD_SHARED_LIBS   OFF CACHE BOOL "" FORCE)
set(BUILD_DOCUMENTATION OFF CACHE BOOL "" FORCE)

FetchContent_Declare(geographiclib
    SOURCE_DIR ${CODE_DIR}/lib/geographiclib
)
FetchContent_MakeAvailable(geographiclib)



# Now we can define our own Libraries by telling the compiler about all the C/C++ files that contain functions that we would like to use.
add_library(spectra ${CODE_DIR}/src/spectra.cpp)
add_library(albedo ${CODE_DIR}/src/albedo.cpp)

add_library(vectorMaths ${CODE_DIR}/src/vectorMaths.cpp)
add_library(ephem ${CODE_DIR}/src/ephem.cpp)
add_library(transforms ${CODE_DIR}/src/transforms/coordTransforms.cpp ${CODE_DIR}/src/transforms/transforms.cpp ${CODE_DIR}/src/transforms/time_conv.cpp)
add_library(integrators ${CODE_DIR}/src/integrators/LSODA.cpp ${CODE_DIR}/src/integrators/RK78.cpp)
add_library(satDynamics ${CODE_DIR}/src/satDynamics.cpp)

# We can also add more include directories to each library specifically rather than globally, for instance
# target_include_directories(stuff
#     ${CODE_DIR}/lib
# )
target_include_directories(integrators PUBLIC
    ${CODE_DIR}/include/integrators
    ${CODE_DIR}/include
    ${CODE_DIR}/lib
)


target_link_libraries(spectra PUBLIC 
    m
)

target_link_libraries(albedo PUBLIC
    m
    netCDF::netcdf
)

target_link_libraries(vectorMaths PUBLIC
    m
)

# If instead we want to add a dependency to an external library that is pre-compiled we add a dependency
# This add dependency guarantees that cspice is built fully before it attempts to build transforms
add_dependencies(transforms cspice_build) 
target_link_libraries(transforms PUBLIC
    m
    vectorMaths
)

add_dependencies(ephem cspice_build)
target_link_libraries(ephem PUBLIC
    m
    GeographicLib::GeographicLib
    ${CSPICE_LIB}
    vectorMaths
    transforms
)

target_link_libraries(satDynamics PUBLIC
    ephem
    vectorMaths
    integrators
    transforms
)




# We can now define the executables that we need to compile, an executable can contain multiple C/C++ files separated by a space.
add_executable(simulation ${CODE_DIR}/src/main.cpp)

# Now we link the executable to all our libraries, it can be private since the executable doesn't need to tell other files about the packages it depends on
target_link_libraries(simulation PRIVATE
    lcm
    spectra
    albedo
    satDynamics
    ephem
    transforms
)