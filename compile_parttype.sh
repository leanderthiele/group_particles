# set to one of GAS, DM, STARS, BH
PartType=$2

g++ -std=c++17 -O3 -ffast-math -funroll-loops \
  -D$PartType \
  -g3 -Wall -Wextra -Wno-unused-parameter -Wno-reorder \
  -I./include -I./include/callback_utils -I./detail \
  -o $1_$PartType $1.cpp \
  -lhdf5 -lhdf5_cpp #-fopenmp
