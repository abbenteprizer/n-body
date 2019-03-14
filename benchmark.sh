#!/bin/bash
g++ -o pnbody pnbody.cpp -fopenmp
g++ -o snbody snbody.cpp

./pnbody 120 10000 1  >  data/bench120.dat
./pnbody 120 10000 2  >>  data/bench120.dat
./pnbody 120 10000 3  >>  data/bench120.dat
./pnbody 120 10000 4  >>  data/bench120.dat
./pnbody 120 10000 5  >>  data/bench120.dat
./pnbody 120 10000 6  >>  data/bench120.dat
./pnbody 120 10000 7  >>  data/bench120.dat
./pnbody 120 10000 8  >>  data/bench120.dat
./pnbody 120 10000 9  >>  data/bench120.dat
./pnbody 120 10000 10  >>  data/bench120.dat
./pnbody 120 10000 11  >>  data/bench120.dat
./pnbody 120 10000 12  >>  data/bench120.dat
./pnbody 120 10000 13  >>  data/bench120.dat
./pnbody 120 10000 14  >>  data/bench120.dat
./pnbody 120 10000 15  >>  data/bench120.dat
./pnbody 120 10000 16  >>  data/bench120.dat

./snbody 120 10000 > data/sbench120.dat

./pnbody 180 10000 1  >  data/bench180.dat
./pnbody 180 10000 2  >>  data/bench180.dat
./pnbody 180 10000 3  >>  data/bench180.dat
./pnbody 180 10000 4  >>  data/bench180.dat
./pnbody 180 10000 5  >>  data/bench180.dat
./pnbody 180 10000 6  >>  data/bench180.dat
./pnbody 180 10000 7  >>  data/bench180.dat
./pnbody 180 10000 8  >>  data/bench180.dat
./pnbody 180 10000 9  >>  data/bench180.dat
./pnbody 180 10000 10  >>  data/bench180.dat
./pnbody 180 10000 11  >>  data/bench180.dat
./pnbody 180 10000 12  >>  data/bench180.dat
./pnbody 180 10000 13  >>  data/bench180.dat
./pnbody 180 10000 14  >>  data/bench180.dat
./pnbody 180 10000 15  >>  data/bench180.dat
./pnbody 180 10000 16  >>  data/bench180.dat

./snbody 180 10000 > data/sbench180.dat


./pnbody 240 10000 1  >  data/bench240.dat
./pnbody 240 10000 2  >>  data/bench240.dat
./pnbody 240 10000 3  >>  data/bench240.dat
./pnbody 240 10000 4  >>  data/bench240.dat
./pnbody 240 10000 5  >>  data/bench240.dat
./pnbody 240 10000 6  >>  data/bench240.dat
./pnbody 240 10000 7  >>  data/bench240.dat
./pnbody 240 10000 8  >>  data/bench240.dat
./pnbody 240 10000 9  >>  data/bench240.dat
./pnbody 240 10000 10  >>  data/bench240.dat
./pnbody 240 10000 11  >>  data/bench240.dat
./pnbody 240 10000 12  >>  data/bench240.dat
./pnbody 240 10000 13  >>  data/bench240.dat
./pnbody 240 10000 14  >>  data/bench240.dat
./pnbody 240 10000 15  >>  data/bench240.dat
./pnbody 240 10000 16  >>  data/bench240.dat

./snbody 240 10000 > data/sbench240.dat
