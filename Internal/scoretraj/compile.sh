#!/bin/bash
g++ -std=c++17 -I $GEMSHOME/gmml/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ scoretraj.cpp -o ../../Bin/scoretraj.exe -lgmml -lpthread 
