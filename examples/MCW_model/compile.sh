#!/bin/sh

clang++ -std=c++11 MCW_model.cpp ../../src/*.cpp ../../src/matrices/*.cpp  -O2 -o ../../bin/MCW_model
