#!/bin/bash

# Script for compiling C++ scripts with the Enzyme automatic differentiation library. Tested on WSL Ubuntu.
# Script must have execute permission to run. This is done by running the following:
#       chmod u+x build.sh 


# First time install of brew (package manager) and Ezyme copiler
# /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
# brew install enzyme

# Compiling Enzyme .cpp file. Explanation found at https://enzyme.mit.edu/getting_started/UsingEnzyme/ 
MAIN_FILE_PATH="main.cpp"
ENZYME_LIBRARY_PATH="/home/linuxbrew/.linuxbrew/Cellar/enzyme/0.0.249_1/lib/LLVMEnzyme-22.so"

clang++ $MAIN_FILE_PATH -std=c++23 -S -emit-llvm -o input.ll -O2 -fno-vectorize -fno-slp-vectorize -fno-unroll-loops
opt input.ll --load-pass-plugin=$ENZYME_LIBRARY_PATH -passes=enzyme -o output.ll -S
opt output.ll -O2 -o output_opt.ll -S
# clang++ output_opt.ll -std=c++23 -DNDEBUG -O3 -o main # prod
clang++ output_opt.ll -std=c++23 -ggdb -pedantic-errors -Wall -Weffc++ -Wextra -Wconversion -Wsign-conversion -Wshadow -Werror -o main # dev

rm input.ll
rm output.ll
rm output_opt.ll

./main
rm ./main
