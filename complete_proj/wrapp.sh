#! /usr/bin/bash

#path to datasets
datasets_path=/mnt/c/Users/Home/OneDrive/uchebnici/Master/3-Semester/SoftwareProj/datasets;

#input alignment file and make a separate folder if there is none
echo "Alignment file: $1";
mkdir ${1%.*}_Analysis #takes only the file name w/out extension


read -p "Run IQTree? (y/n)" iqtr;
case ${iqtr:0:1} in
    y|Y )  #run iqtree to get the ML tree, implement later
        echo yes
    ;;
    * )   #only run the tests
        ./a -F $datasets_path/$1/alignment.nex > ${1%.*}_Analysis/${1%.*}_full.txt
    ;;
esac









#ggf model input for iqtree
