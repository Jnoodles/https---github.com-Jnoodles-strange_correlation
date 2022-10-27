#!/bin/bash
g++ same_kpom_kpop_deduction.C -o kpop  `root-config --cflags --glibs`
echo "kpom-kpop"
sleep 5m
./kpop 1 11 1
