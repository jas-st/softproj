#!/usr/bin/env bash

#exits with an error
set -e
start=$SECONDS;
iqtr=false;
s=false;

#set flag options & variables
while getopts ":a:t:m:s:I" flag
do
    case "${flag}" in
        a) alignment_file=${OPTARG};;
        t) treefile=${OPTARG};;
        m) model=${OPTARG};;
        s) s=${OPTARG};;
        I) iqtr=true;;
      ##ERROR: Please provide an alignment file
    esac
done

#check if alignment file provided
if [ -z "$alignment_file" ]; then
    echo "You did not specify an alignment file!";
    echo "USAGE: bash analysis_biological.sh -a alignment_file.nex"
    echo "exiting..."
    exit
fi

#check if alignment file exists
if [ ! -f "$alignment_file" ]; then
    echo "Alignment file $alignment_file does not exist."
    echo "exiting ...";
    exit
fi

#set additional variables and default options
alignment_name=${alignment_file%.*};
alignment_name=${alignment_name##*/};


#print messages
echo "Alignment file: " $alignment_file;
echo "Saturation tests: " $s;
echo "Results in directory: results_"$alignment_name;

mkdir -p results_$alignment_name;
rm -f results_$alignment_name/*.pdf;
rm -f results_$alignment_name/*.csv;
rm -f results_$alignment_name/*.tree;
rm -f results_$alignment_name/*.log;

if [ -z "$treefile" ]; then
    if [ "$iqtr" = false ] ; then
    	echo -e "No tree file/IQ-TREE option provided. Only computing the test statistics."
    else
    	mkdir -p results_$alignment_name/IQTree_Results
    	cp $alignment_file results_$alignment_name/IQTree_Results/;
    	cd results_$alignment_name/IQTree_Results;
	
    	if [ -z "$model" ] ; then
    	    echo "Using ML-Tree computed with IQ-TREE and ModelFinder."
    	    echo "-------------------------------"
    	    echo "...Running IQ-TREE..."
    	    #touch ${alignment_file##*/}.treefile
    	    iqtree2 -s "${alignment_file##*/}" --redo-tree >> tree_inference.log;
    	else
    	    echo "Using ML-Tree computed with IQ-TREE and the "${model}" model."
    	    echo "-------------------------------"
    	    echo "...Running IQ-TREE..."
    	    iqtree2 -s "${alignment_file##*/}" -m "${model}" --redo-tree >> tree_inference.log;
    	    #touch ${alignment_file##*/}.treefile
    	fi
    	treefile=results_$alignment_name/IQTree_Results/${alignment_file##*/}.treefile;
    	echo "Finished!"
    	cd ../..;
    fi
else
    if [ ! -f "$treefile" ]; then
       echo "Tree file $treefile does not exist."
       echo "exiting ...";
       exit
    fi
    echo "Using tree: " $treefile;
fi

#check operation system for mac (needs different c++ file)
cfile=all_tests.out

if [[ "$OSTYPE" =~ ^darwin ]]; then
    cfile=all_tests_mac.out
fi

#C++ script - Test Statistics computation
echo "-------------------------------------"
echo "...Calculating test statistics..."
./bin/$cfile -F "$alignment_file" -s $s> results_$alignment_name/results_raw_${alignment_name}.csv;
#touch ./results_$alignment_name/${alignment_name}_tests.csv;
echo "Computation done!"

if [ -n "$treefile" ]; then
	#R Script - Analysis & Visualisation
	echo "-------------------------------------"
	echo "...Running R Script..."
	Rscript ./bin/analysis_visualisation_biological.R $treefile results_$alignment_name/results_raw_${alignment_name}.csv $s >> results_$alignment_name/r_log.log 2> >(tee results_$alignment_name/r_log.log >&2);
	#touch ./results_$alignment_name/${alignment_name}_tests.pdf;
	echo "Computation done!"
fi


rm -f Rplots.pdf;


echo "Duration in seconds:" $(( SECONDS - start ));
