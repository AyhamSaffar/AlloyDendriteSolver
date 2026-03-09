cmake_minimum_required(VERSION 3.20)

find_package(Enzyme REQUIRED CONFIG) # must install enzyme & lld through Brew / Spack
set(CMAKE_CXX_COMPILER "${Enzyme_LLVM_BINARY_DIR}/bin/clang++") # enzyme only on the compiler it was installed on
set(CMAKE_CXX_STANDARD 20) # requires at least C++ 20

project(AlloySolidificationSolver)
add_executable(AlloySolidificationSolver src/main.cpp)
target_include_directories(AlloySolidificationSolver PUBLIC include) # headers in include directory
target_link_libraries(AlloySolidificationSolver PUBLIC LLDEnzymeFlags)
