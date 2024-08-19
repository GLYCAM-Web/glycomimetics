#!/bin/bash
g++ -std=c++17 -I $GEMSHOME/gmml/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ -g -fvar-tracking scoretraj.cpp -o scoretraj.exe -lgmml -lpthread 
#g++ -std=c++17 -I $GEMSHOME/gmml/ -L$GEMSHOME/gmml/bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ -g -fvar-tracking main.cpp -lgmml -pthread -o ../main.exe
#g++ -std=c++11 -I$GEMSHOME/gmml/includes -L$GEMSHOME/gmml/lib -L/usr/include/c++/9/boost_1_66_0/stage/lib make_ligand_pdb2entropy_metadata.cpp -o make_ligand_pdb2entropy_metadata.out -lgmml -lpthread -lboost_filesystem -lboost_system 
#g++ -std=c++11 -I$GEMSHOME/gmml/includes -L$GEMSHOME/gmml/lib -L/usr/include/c++/9/boost_1_66_0/stage/lib rename_ligand.cpp -o rename_ligand.out -lgmml -lpthread -lboost_filesystem -lboost_system 
