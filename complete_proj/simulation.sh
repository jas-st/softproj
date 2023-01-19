#!/usr/bin/env bash

treefile=$1;
treename=${treefile%.*};
seq_len=${3:-500}
matrix=${2:-'JC'}
k=${4:-1}

echo "This is tree name: $treename";

echo '' > results_raw_$treename.tsv

mkdir -p alignment_alisim;
cp $treefile alignment_alisim/;

for value in $(seq $k);
do
	echo "Simualtion: " $value;	
	cd alignment_alisim;
	iqtree2 --alisim $treename-alignment -m $matrix -t $1 --length $seq_len > /dev/null;
	cd ..;
	./all_tests.out -F alignment_alisim/$treename-alignment.phy >> results_raw_$treename.tsv;
done

mkdir -p plots;
Rscript analysis.R $treefile $seq_len results_raw_$treename.tsv results_$treename.csv $k > /dev/null; 


echo "Tree file: " $treefile;
echo "Length of simulated sequences: " $seq_len;
echo "Model of substitution: "$matrix;
echo "Number of simulations: "$k;

