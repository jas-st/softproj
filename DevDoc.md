# Developer Level Documentation

## Bash 

### `analysis_simulation.sh`

##### Syntax 
```
bash analysis_simulation.sh -t <tree> [OPTIONS]
```

##### Required parameters
- `-t` - tree file in Newick format

##### Options
- `-m` - substitution model for sequence simulation, by default JC 
- `-n` - sequence length, by default 500
- `-k` - number of simulations, by default 1
- `-s` - apply staturation tests true or false, by default false

##### Usage 
```
bash analysis_simulation.sh -t tree.nwk -m JC -n 1000 -k 100  -s true
```

##### Structure
- Call Alisim
```
iqtree2 --alisim $treename-alignment -m $matrix -t ${t} --length $seq_len --num-alignments $k -seed 123 >> $dir/results_$treename/
```
- Call C++ Script
- Call R Script

### `analysis_biological.sh`

##### Synthax
```
bash analysis_biological.sh -a <alignment> [OPTIONS]
```

##### Required parameters
- `-a` - multiple sequence alignment in Phylip or NEX format

##### Options
- `-t`
- `-iqtr`
- `-m`
- `-s`

##### Structure
- Call C++ Script
- Call IQ-TREE 
- Call R Script

## C++

### *Sequence*
Class for storing sequences.   
Location: `bin\lib\Sequence.h` 

##### Constructors 
```
Sequence(std::string seq_id, std::string seq_str, Eigen::Vector4d seq_freq)
```

##### Parameters
- `id` - sequence id/name/etc
- `seq` - nucleotide sequnce
- `nucl_freqs` - vector

##### Methods


### *Alignment*
Class for storing sequences  

##### Syntax
Location: `bin\lib\Sequence.h`  
```
#include<'Sequence.h'>

Sequence(std::string seq_id, std::string seq_str, Eigen::Vector4d seq_freq)
```

##### Parameters
- `id` - sequence id/name/etc
- `seq` - nucleotide sequnce
- `nucl_freqs` - vector

##### Methods

### *seqs_read*

Function for reading in sequences in an Alignment object.

#### Usage
```
Alignment seqs_read(std::string file_name, std::string extension)
```

#### Parameters


## R

### *head_success*
Function for plotting a heatmap to visualize the pairwise p-values and decision results. 

##### Usage 
```
heat_sucess(seq_pair, test)
```

##### Parameters
- `seq_pair` - data frame column or character vector that contains the sequence pairs
- `test_pv`- data frame column or numeric vector that contains the p-values

### *edges_rej*
Function for mapping the pairwise test results on the tree. Returns a dataframe with edge id (from the ggtree object) in the first column and frequency of rejections in the second.

##### Usage 
```
edjes_rejected_freq(tree, seq_pair, test_pv)
```

##### Parameters
- `tree`- ggtree object
- `seq_pairs` - data frame column or character vector that contains the sequence pairs
- `test_pv`- data frame column or numeric vector that contains the p-values

### *coloured_tree*
Function for plotting a coloured tree. Calculates edges_rej data frame and uses the values as colour and label in the treeplot.

##### Usage 
```
coloured_tree(tree, seq_pair, test_pv)
```

##### Parameters
- `tree`- ggtree object
- `seq_pairs` - data frame column or character vector that contains the sequence pairs
- `test_pv`- data frame column or numeric vector that contains the p-values

