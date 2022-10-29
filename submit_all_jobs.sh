#!/bin/bash

jobDir=$PWD
a=1
for i in {1..171}
do
    cd ${jobDir}
    if [ ! -d "$i" ]; then
        mkdir $i
    fi
    cp slurm_1.sh slurm_2.sh $i/
    cp command_1.sh command_2.sh $i/
    cp same_kmom_kmop_deduction.C same_kpom_kpop_deduction.C $i/
    cd ${jobDir}/$i
    sed -i "s/name=ko/name=ko${i}1/g" slurm_1.sh
    sed -i "s/name=ko/name=ko${i}2/g" slurm_2.sh
    sed -i "s/.\/kmom/.\/kmom $a $[$a+10] $i/g" command_1.sh
    sed -i "s/.\/kpop/.\/kpop $a $[$a+10] $i/g" command_2.sh
    sbatch slurm_1.sh
    sbatch slurm_2.sh
    echo $a
    let a=a+10

done
