#!/bin/bash
jobDir=$PWD
g++ same_kmom_kmop_deduction.C -o kmom  `root-config --cflags --glibs`
echo "kmom-kmop"
./kmom 1 11 1
