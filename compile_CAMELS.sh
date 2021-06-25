# set to one of ILLUSTRIS, SIMBA
SimType=$2

g++ -std=c++17 -O3 -ffast-math -funroll-loops \
  -DFOR_$SimType \
  -g3 -Wall -Wextra -Wno-unused-parameter -Wno-reorder \
  -I./include -I./include/callback_utils -I./detail \
  -o $1_$SimType $1.cpp \
  -lhdf5 -lhdf5_cpp -fopenmp
