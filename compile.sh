module load hdf5/gcc/1.10.0
module load rh/devtoolset/8

g++ -std=c++17 -O3 -ffast-math -g3 -Wall -Wextra -Wno-unused-parameter -Wno-reorder -o $1 $1.cpp -lhdf5 -lhdf5_cpp
