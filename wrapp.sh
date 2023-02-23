#! /usr/bin/bash
set -o errexit

#input alignment file
file=$1;

#extract name
name=${file%.*};
name=${name##*/}

#extract extension
ext=${file##*.}

#also possible to input tree name
tree=$2;

#verify if file exists and has the correct extension
if ! test -f "$file";
then
    echo "Error: File doesn't exist."
    exit 1
elif [ "$ext" != "nex" ] && [ "$ext" != "phy" ];
then
    echo "Error: Please provide a .nex or a .phylip alignment file."
    exit 1
fi

#create a directory where the analysis file are going to be stored
echo -e "Alignment file: $file \n";
#maybe ask if the user wants to override it??
mkdir -p results_$name;
cp $file results_$name;

#run the tests regardless if we have the tree
#echo "...Calculating test statistics..."
#./a -F "$file" > ${name}_tests.csv;
#mv ${name}_tests.csv results_$name;
#echo "Done! Created CSV file in directory results_$name."


#check if the user has provided a tree file
if [ -z "$tree" ] ;
then
    #if tree file name is empty ask for IQtree
    read -p "Run IQTree? (y/n)" iqtr;
    case ${iqtr:0:1} in
    y|Y )  #run iqtree to get the ML tree
        echo "...Running IQ-Tree..."
        ../iqtree2.exe -s ../$file --redo-tree -quiet;
        #ggf model input for iqtree
        ;;

    n|N )
        #ask for a tree file if no iqtree
        read -p "Please provide a tree file:" tree;
    esac
fi

if [ -z "$tree" ] ;
then
    echo "Error : No tree file provided. Cannot run analysis."
fi

#if tree run R script
    #R script analysis
