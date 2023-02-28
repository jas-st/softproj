# Developer Level Documentation

## Bash 

### `analysis_simulation.sh`
Bash script that acts like a wrapper around the different functions of this workflow. 

##### Usage 
```
bash analysis\_simulation.sh -t <treename>.nwk -m JC -n 1000 -k 100  -s true
```

##### Parameters
- `t` - tree file in Newick format `<treename>.nwk`
- `m` - substitution model for sequence simulation `<m>`, by default JC 
- `n` - sequece length `<n>`, by default 500
- `k` - number of simulations `<k>`, by default 1
- `s` - apply staturation tests true or false `<s>`, by default false

##### Structure
1. read in the command line arguments
2. create directory results_<treename>, where all results will be saved
3. create a new directory `alignment_alisim`
4. call the IQ-Tree command (save command line output to `simulation.log`): 

```
  iqtree2 –alisim alignment -m <model> -t <treename>.nwk –length <n> –num-alignments <k> 
```
5. loop trough all alignment simulation files in `alignment_alisim`
6. apply C++ script by using command:
```
  ./all_tests.out -F <treename>.nwk -s <s>
```
7. save raw test results in results_<treename>/results_raw_<treename>.csv
8. call R script on the raw test results:
```
  analysis.R <treename>.nwk <n> results_<treename>/results_raw_<treename>.csv <k> <s>
```
9. delete directory with simulated alignments `alignment_alisim`
  
  
### `analysis_biological.sh`

## C++

### *Sequence*
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


## R file: analysis_visualisation_simulation.R / analysis_visualisation_biological.R
Used for calculating p-values and decisions of test statistics and plotting the results with 4 different visualsations. 
  
#### Usage
```
analysis.R <treefile> <n> <raw_test_statistcs> <k> <s>
```  
#### Parameters
- `treefile` - treefile either in Newick format or 
- `n`- sequence length
- `raw_test_statistcs` - path to .csv file with raw test statistics
- `k` - number of simulations
- `s` - true/false saturation test

Includes functions:

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
  
