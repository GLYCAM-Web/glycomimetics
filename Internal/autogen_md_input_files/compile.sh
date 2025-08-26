#!/bin/bash
g++ -std=c++17 -I $GEMSHOME/gmml/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ main_gbsa.cpp -o ../../Bin/make_gbsa_input.exe -lgmml -lpthread 
g++ -std=c++17 -I $GEMSHOME/gmml/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ main_pbsa.cpp -o ../../Bin/make_pbsa_input.exe -lgmml -lpthread 
