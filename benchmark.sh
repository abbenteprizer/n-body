#!/bin/bash
g++ -o pnbody pnbody.cpp -fopenmp
gnuplot imageScriptParallel.p
