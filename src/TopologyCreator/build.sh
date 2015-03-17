#!/bin/sh

g++ -I/usr/local/cuda/include main.cpp pdbio.cpp psfio.cpp topio.cpp aatocg.cpp ../Util/parameters.cpp ../Util/wrapper.cpp ../IO/configreader.cpp -o sop-top2
