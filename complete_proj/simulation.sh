#!/usr/bin/env bash

treefile=$1;
treename=${treefile%.*};
seq_len=${3:-500}
matrix=${2:-'JC'}
k=${4:-1}

echo "This is tree name: $treename";

mkdir -p alignment_alisim;
cp $treefile alignment_alisim/;
mkdir -p results_$treename;
echo '' > results_$treename/results_raw_$treename.tsv;

for value in $(seq $k);
do
	echo "Simualtion: " $value;	
	cd alignment_alisim;
	iqtree2 --alisim $treename-alignment -m $matrix -t $1 --length $seq_len > /dev/null;
	cd ..;
	./all_tests_new.out -F alignment_alisim/$treename-alignment.phy >> results_$treename/results_raw_$treename.tsv;
done

Rscript analysis.R $treefile $seq_len results_$treename/results_raw_$treename.tsv results_$treename/results_$treename.csv $k; 


echo "Tree file: " $treefile;
echo "Length of simulated sequences: " $seq_len;
echo "Model of substitution: "$matrix;
echo "Number of simulations: "$k;
echo "Results in directory: results_"$treename; 
