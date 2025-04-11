cmake_minimum_required(VERSION 3.15)
project(BTR_Field_3D VERSION 0.1)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Pybind11 integration
find_package(Python REQUIRED COMPONENTS Interpreter Development)
add_subdirectory(extern/pybind11)

# --------------------------------------------------
# Library: Interpolation
# --------------------------------------------------
add_library(Interpolation
  src/interpolation/interpolation.cpp
  src/interpolation/interpolation_constants.cpp
)

target_include_directories(Interpolation
  PUBLIC 
    include
  PRIVATE
    src/interpolation
)

# --------------------------------------------------
# Library: MF
# --------------------------------------------------
add_library(MF
  src/mf/read_field.cpp
  src/mf/get_field.cpp
  src/mf/interpolate_field.cpp
)

target_include_directories(MF
  PUBLIC 
    include
  PRIVATE
    src/mf
)

target_link_libraries(MF
  PUBLIC 
    Interpolation
)

# --------------------------------------------------
# Python Bindings
# --------------------------------------------------
pybind11_add_module(mf_3d
  python/mf_python.cpp
)

target_link_libraries(mf_3d
  PRIVATE
    MF
    Interpolation
)

# --------------------------------------------------
# Main Executable
# --------------------------------------------------
add_executable(BTR_Field_3D
  src/main.cpp
)

target_link_libraries(BTR_Field_3D
  PRIVATE
    MF
    Interpolation
)
