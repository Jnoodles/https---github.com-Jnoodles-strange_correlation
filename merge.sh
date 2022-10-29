#!/bin/bash

#jobdir=$PWD
#echo ${jobdir}

singularity shell ~/pyroot.sif
jobdir=$(dirname "$PWD")
echo $jobdir
pair1=0
pair2=0
pair3=0
pair4=0
echo -e "\t K-O- \t K-O+ \t K+O- \t K+O+"
for i in {1..170}
do

    if [ ! -d "${jobdir}/$i" ]; then
        continue
    fi
    cd ${jobdir}/$i/

    eve1=`cat output-*1 | grep "same pair 1 events: " | cut -c 21-40`
    eve2=`cat output-*1 | grep "same pair 2 events: " | cut -c 21-40`
    eve3=`cat output-*2 | grep "same pair 1 events: " | cut -c 21-40`
    eve4=`cat output-*2 | grep "same pair 2 events: " | cut -c 21-40`
    echo -e "$i\t $eve1 \t $eve2 \t $eve3 \t $eve4"
    let pair1+=eve1
    let pair2+=eve2
    let pair3+=eve3
    let pair4+=eve4
done
echo -e "total \t $pair1  $pair2  $pair3  $pair4"
cd ${jobdir}/all
#plot

hadd kmom_kmop_dong_default.root ../{1..171}/deduction_kmom_kmop_*.root
hadd kpom_kpop_dong_default.root ../{1..171}/deduction_kpom_kpop_*.root
root "merge_ko.C+($pair1, $pair2,$pair3, $pair4)"
#echo 'exit'

