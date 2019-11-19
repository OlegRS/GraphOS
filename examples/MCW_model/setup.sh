#!/bin/sh

[ -d "../../data" ] && ([ -d "../../data/MCW_example/" ] || mkdir ../../data/MCW_example/) || (mkdir  ../../data && mkdir ../../data/MCW_example/)

[ -d ../../bin ] || mkdir  ../../bin

clang++ -std=c++11 MCW_model.cpp ../../src/*.cpp ../../src/matrices/*.cpp  -O2 -o ../../bin/MCW_model
