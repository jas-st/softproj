#!/usr/bin/env bash

set -e
start=$SECONDS

while getopts t:n:m:k:s: flag
do
    case "${flag}" in
        t) t=${OPTARG};;
        n) n=${OPTARG};;
        m) m=${OPTARG};;
        k) k=${OPTARG};;
        s) s=${OPTARG};;
        ##ERROR: Please provide an tree file
    esac
done

#set additional variables and default options
treefile=${t};
treename=${treefile%.*};
seq_len=${n:-500}
matrix=${m:-'JC'}
k=${k:-1}
dir=$(pwd)

mkdir -p results_$treename;
rm -f results_$treename/* || true;

#logfile
touch results_$treename/simulation.log;
logfile=$dir/results_$treename/simulation.log;

echo "Tree file: " ${t} | tee -a $logfile;
echo "Length of simulated sequences: " $seq_len | tee -a $logfile;
echo "Model of substitution: "$matrix | tee -a $logfile;
echo "Number of simulations: "$k | tee -a $logfile;
echo "Results in directory: results_"$treename | tee -a $logfile;
echo "This is tree name: $treename" | tee -a $logfile;
echo "Saturation test:" ${s} | tee -a $logfile;


#check if treefile exists
if [ ! -f "$treefile" ]; then
    echo "$treefile does not exist."
fi

nwktree=$(<"$treefile");
if [[ $nwktree =~ "aaa" ]]; then
   echo "Replace aaa in" $treefile "with" $m "!" | tee -a $logfile;
   nwktree=${nwktree/aaa/$m};
   echo "New Newick string: "$nwktree | tee -a $logfile;
fi

echo $nwktree >> results_$treename/$treefile;
mkdir -p alignment_alisim;
rm -f alignment_alisim/* || true;
cp results_$treename/$treefile alignment_alisim/;


echo "Simulating alignment with IQ-TREE alisim."
echo "-------------------------------"
echo "...Running IQ-TREE alisim..."
echo '' > results_$treename/results_raw_$treename.csv;
cd alignment_alisim;
iqtree2 --alisim $treename-alignment -m $matrix -t ${t} --length $seq_len --num-alignments $k -seed 123  | tee -a $logfile;
cd ..;

echo "Simulations Finished!";

#C++ script - Test Statistics computation
echo "-------------------------------------"
echo "...Calculating test statistics..."


#check operation system for mac (needs different c++ file)
cfile=all_tests.out

if [[ "$OSTYPE" =~ ^darwin ]]; then
    cfile=all_tests_mac.out
fi

for f in ./alignment_alisim/$treename-alignment*.phy;
do
	echo "Test on simulation file: "$f;	
	./bin/$cfile -F $f -s ${s} >> results_$treename/results_raw_$treename.csv;
done
echo "Computation Finished!"

#R Script - Analysis & Visualisation
echo "-------------------------------------"
echo "...Running R Script..."
Rscript bin/analysis.R $treefile $seq_len results_$treename/results_raw_$treename.csv results_$treename/results_$treename.csv $k ${s} | tee -a $logfile;
rm -r alignment_alisim;
mv Rplots.pdf results_$treename/venn_diag.pdf
echo "Analysis and Visualisation Finished!"

echo "Duration in seconds:" $(( SECONDS - start )) | tee -a $logfile;
