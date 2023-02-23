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
    esac
done

treefile=${t};
treename=${treefile%.*};
seq_len=${n:-500}
matrix=${m:-'JC'}
k=${k:-1}
dir=$(pwd)

echo "Tree file: " ${t};
echo "Length of simulated sequences: " $seq_len;
echo "Model of substitution: "$matrix;
echo "Number of simulations: "$k;
echo "Results in directory: results_"$treename;
echo "This is tree name: $treename";
echo "Saturation test:" ${s};

echo $dir;
mkdir -p alignment_alisim;
rm -f alignment_alisim/* || true;
cp $treefile alignment_alisim/;
mkdir -p results_$treename;
rm -f results_$treename/* || true;

#logfile
touch results_$treename/simulation.log;

echo '' > results_$treename/results_raw_$treename.csv;
cd alignment_alisim;
iqtree2 --alisim $treename-alignment -m $matrix -t ${t} --length $seq_len --num-alignments $k -seed 123 >> $dir/results_$treename/simulation.log;
cd ..;


echo "Simulations Done!";


for f in ./alignment_alisim/$treename-alignment_*.phy;
do
	echo "Test on simulation file: "$f;	
	./bin/all_tests.out -F $f -s ${s} >> results_$treename/results_raw_$treename.csv;
done

Rscript bin/analysis.R $treefile $seq_len results_$treename/results_raw_$treename.csv results_$treename/results_$treename.csv $k ${s} --quiet >> $dir/results_$treename//simulation.log; 
rm -r alignment_alisim;

echo "Duration in seconds:" $(( SECONDS - start )); 
