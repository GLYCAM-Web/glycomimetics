#!/bin/bash
g++ -std=c++17 -I $GEMSHOME/gmml/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ main.cpp -o ../../Bin/make_resp_in.exe -lgmml -lpthread 
