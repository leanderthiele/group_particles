g++ -std=c++17 -O3 -ffast-math -funroll-loops \
  -g3 -Wall -Wextra -Wno-unused-parameter -Wno-reorder \
  -I./include -I./include/callback_utils -I./detail \
  -o $1 $1.cpp \
  -lhdf5 -lhdf5_cpp -fopenmp
