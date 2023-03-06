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

#set additional variables and default options
alignment_name=${alignment_file%.*};
alignment_name=${alignment_name##*/};


#print messages
echo "Alignment file: " $alignment_file;
echo "Saturation tests: " $s;
echo "Results in directory: results_"$alignment_name;

mkdir -p results_$alignment_name;
rm -f results_$alignment_name/*.pdf || true;
rm -f results_$alignment_name/*.csv || true;
rm -f results_$alignment_name/*.tree || true;

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
    echo "Using tree: " $treefile;
fi

#C++ script - Test Statistics computation
echo "-------------------------------------"
echo "...Calculating test statistics..."
./bin/all_tests.out -F "$alignment_file" -s $s> results_$alignment_name/results_raw_${alignment_name}.csv;
#touch ./results_$alignment_name/${alignment_name}_tests.csv;
echo "Computation done!"

if [ -n "$treefile" ]; then
	#R Script - Analysis & Visualisation
	echo "-------------------------------------"
	echo "...Running R Script..."
	Rscript ./bin/analysis_visualisation_biological.R $treefile results_$alignment_name/results_raw_${alignment_name}.csv $s >> r_log.log;
	#touch ./results_$alignment_name/${alignment_name}_tests.pdf;
	echo "Computation done!"
fi

echo "Duration in seconds:" $(( SECONDS - start ));
