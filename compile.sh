module load hdf5/gcc/1.10.0
module load rh/devtoolset/8

g++ -std=c++17 -O3 -g3 -Wall -Wextra -Wno-unused-parameter Y_Delta.cpp -lhdf5 -lhdf5_cpp