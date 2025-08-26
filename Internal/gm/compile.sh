#!/bin/bash
g++ -std=c++17 -I $GEMSHOME/gmml/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ main.cpp -o ../../Bin/gm.exe -lgmml -lpthread 
#g++ -std=c++17 -I $GEMSHOME/gmml/ -static -static-libgcc -static-libstdc++ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ main.cpp -o ../../Bin/gm_static.exe -lgmml -lpthread -lc -lm
