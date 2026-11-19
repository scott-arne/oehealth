cmake_minimum_required(VERSION 3.16)

if(POLICY CMP0144)
    cmake_policy(SET CMP0144 NEW)
endif()

project(oehealth VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Options
option(OEHEALTH_BUILD_TESTS "Build tests" ON)
option(OEHEALTH_BUILD_PYTHON "Build Python bindings" ON)
option(OEHEALTH_UNIVERSAL2 "Build universal2 binary for macOS" OFF)
option(OEHEALTH_USE_STABLE_ABI "Use Python stable ABI" ON)

if(APPLE AND OEHEALTH_UNIVERSAL2)
    set(CMAKE_OSX_ARCHITECTURES "x86_64;arm64" CACHE STRING "Build architectures" FORCE)
elseif(APPLE)
    unset(CMAKE_OSX_ARCHITECTURES CACHE)
endif()

if(APPLE AND NOT CMAKE_OSX_DEPLOYMENT_TARGET)
    set(CMAKE_OSX_DEPLOYMENT_TARGET "11.0" CACHE STRING "Minimum macOS version")
endif()

# Dependencies
include(FetchContent)

FetchContent_Declare(
    cmake_openeye
    GIT_REPOSITORY https://github.com/scott-arne/cmake-openeye.git
    GIT_TAG v1.2.2
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(cmake_openeye)
list(APPEND CMAKE_MODULE_PATH "${cmake_openeye_SOURCE_DIR}")

if(OEHEALTH_BUILD_PYTHON)
    find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
    find_package(SWIG 4.0)
endif()

if(OPENEYE_USE_SHARED)
    find_package(OpenEyePython REQUIRED)
endif()

find_package(OpenEye REQUIRED)

# nlohmann/json is header-only and used only by src/json.cpp.
FetchContent_Declare(
    nlohmann_json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG v3.11.3
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(nlohmann_json)

# Library
set(OEHEALTH_SOURCES
    src/assessment.cpp
    src/assessments/bond_angles.cpp
    src/assessments/bond_lengths.cpp
    src/assessments/clashes.cpp
    src/assessments/geometry_outliers.cpp
    src/assessments/ramachandran.cpp
    src/assessments/valence.cpp
    src/common.cpp
    src/context.cpp
    src/data/engh_huber_2001.cpp
    src/geometry/measure.cpp
    src/geometry/targets.cpp
    src/json.cpp
    src/result.cpp
    src/suite.cpp
)

add_library(oehealth ${OEHEALTH_SOURCES})

set_target_properties(oehealth PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(oehealth
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# OERamachandranAnalysis lives in OEBio, which depends on OEGrid. The
# BUILD_INTERFACE wrapper keeps the header-only json target out of the
# installed export set.
target_link_libraries(oehealth
    PUBLIC
        OpenEye::OEBio
        OpenEye::OEChem
    PRIVATE
        $<BUILD_INTERFACE:nlohmann_json::nlohmann_json>
)

# OpenEye pulls in windows.h, whose wingdi.h defines an ERROR macro that breaks
# Severity::ERROR and Status::ERROR; NOGDI suppresses it. PUBLIC so the SWIG
# module and the tests inherit it.
target_compile_definitions(oehealth PUBLIC $<$<PLATFORM_ID:Windows>:NOGDI>)

# Tests
if(OEHEALTH_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Python bindings
if(OEHEALTH_BUILD_PYTHON AND SWIG_FOUND)
    add_subdirectory(swig)
endif()

# Install
install(TARGETS oehealth
    EXPORT oehealth-targets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(DIRECTORY include/oehealth DESTINATION include)

install(EXPORT oehealth-targets
    FILE OEHEALTHTargets.cmake
    NAMESPACE OEHEALTH::
    DESTINATION lib/cmake/oehealth
)
