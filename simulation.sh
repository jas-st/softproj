#!/usr/bin/env bash

echo "Use Tree file: $1";

mkdir -p alignment_alisim;
cp $1 alignment_alisim/;
cd alignment_alisim;
iqtree2 --alisim alignment -m JC -t $1 --length 500 > nul;
cd ..;
./a.out -F alignment_alisim/alignment.phy;
