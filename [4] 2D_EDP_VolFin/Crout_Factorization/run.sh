#!/bin/bash

g++ 2D_EDP_CF.cpp -lm -o 2D_EDP_CF && ./2D_EDP_CF

source bin/activate

python3 analise.py