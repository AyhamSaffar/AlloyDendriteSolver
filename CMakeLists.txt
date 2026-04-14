cmake_minimum_required(VERSION 3.20)

find_package(Enzyme REQUIRED CONFIG) # must install enzyme & lld through Brew / Spack
set(CMAKE_CXX_COMPILER "${Enzyme_LLVM_BINARY_DIR}/bin/clang++") # enzyme only works on the compiler it was installed on
set(CMAKE_CXX_STANDARD 20) # required for #include <numbers>

project(AlloySolidificationSolver)

# saves all .cpp files in src directory to SRC_FILES variable. file names must not have spaces! 
# CONFIGURE_DEPENDS reruns CMAKE when new .cpp files are found in src. Might not be supported on all platforms
file(GLOB SCRIPT_FILE_PATHS CONFIGURE_DEPENDS scripts/*.cpp)

foreach(SCRIPT_FILE_PATH ${SCRIPT_FILE_PATHS})
    get_filename_component(SCRIPT_FILE_NAME ${SCRIPT_FILE_PATH} NAME_WE) #remove parent directories & file extension
    add_executable(${SCRIPT_FILE_NAME} ${SCRIPT_FILE_PATH})
    target_include_directories(${SCRIPT_FILE_NAME} PRIVATE include) # headers in include directory
    target_link_libraries(${SCRIPT_FILE_NAME} PRIVATE LLDEnzymeFlags)
    target_compile_definitions(${SCRIPT_FILE_NAME} PRIVATE DATA_PATH="${CMAKE_CURRENT_SOURCE_DIR}/data/")
endforeach(SCRIPT_FILE_PATH ${SCRIPT_FILE_PATHS})

file(GLOB TEST_FILE_PATHS CONFIGURE_DEPENDS tests/*.cpp)

add_subdirectory(libs/Catch2)
add_executable(tests ${TEST_FILE_PATHS})
target_include_directories(tests PRIVATE include)
target_link_libraries(tests PRIVATE Catch2::Catch2WithMain LLDEnzymeFlags)
