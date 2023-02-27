# User Level Documentation
### _Application of a test for reversibility of Markov chains in molecular evolution_

Pipeline for the application of the tests Bowker, Stuart, Internal Symmetry and Quasi-symmetry on simulated or biological data. The software is executable on the command line. It comprises two bash scripts, which serve as a link between a C++ program and an R Script for the analysis of the provided data.
The program uses IQ-Tree to simulate an alignment (AliSim) provided a tree topology, or to find the maximum likelihood tree (ModelFinder) provided an alignment file, which it then uses for the analysis.

## Background
We assume a multiple sequence alignment (MSA) or a phylogenetic tree follows a Markovian model of sequence evolution. A common assumption
when performing phylogenetic analysis based on is that the model is stationary and reversible. This pipline applies 4 different test:

- Bowker Test for Symmetry: stationarity + reversibility if Null hypothesis is kept
- Stuart Test for Symmetry: stationarity if Null hypothesis is kept
- Test for Quasi-Symmetry: reversibility if Null hypothesis is kept
- Test for Internal-Symmetry

Each test is applied pairwise on the MSA. There is not yet a decision rule for the model on a bigger tree, but we aim to help with the pairwise analysis and visualisations. To avoid including non-informative pairs in the decision process, saturation tests can be applied on the MSA. These are applied pairwise again:

- Cassius' Test for Saturation 1: saturation if Null hypothesis is rejected (uses relative nucleotide frequency of the whole MSA)
- Cassius' Test for Saturation 2: saturation if Null hypothesis is rejected (uses relative nucleotide frequency of the sequence pair)
- Chi-square Test for Saturation: saturation if Null hypothesis is rejected

## Structure

softproj\
│   `analysis_simulation.sh` - Bash script for the analysis of simulated data\
│   `analysis_biological.sh` - Bash script for the analysis of real data\
└───bin\
│--------   │   `analysis.R` - the R script\
│--------   │   `all_tests.cpp` - the C++ script\
│--------   │   `all_tests.out` - the compiled C++ script\
│--------   └───lib - all libraries and header files needed for the C++ script\
│----------------       │   headerfile1.h\
│----------------       │   headerfile2.h\
│----------------       │   headerfile3.h\
└───test_input

## Installation & Dependencies

The software requires **R** and **IQ-Tree** to run. 

### R Script - built under R 4.2.2

R must be installed as well as the following packages:
- scales
- tidyverse
- ape
- Biostrings
- ggtree
- ggpubr
- ggvenn
- tidytree
- treeio

### IQ-Tree

The software uses the AliSim extension, as well as ModelFinder, contained in the IQ-TREE software. The user has to have IQ-TREE with AliSim installed and be able to run IQ-TREE by simply entering `iqtree`.(or include iqtree executable in the software??? or provide an option of specifiying the iqtree path)

There is otherwise no installation needed, simply run the desired bash script in the command line, providing the needed input files.

## Usage
To start, open a terminal (cmd/Powershell in Windows) and navigate to the path of the program.

#### `analysis_simulation.sh`
The analysis of simulated data requires a tree file and the parameters for the simulation.
```sh
bash analysis_simulation.sh -t treefile.tree -m JC -n 1000 -k 100 -s true
```
- -t - the tree file in standard Newick format
- -m - specifies a substitution model to use for the simulation (default: JC), all possible options can be seen in [substitution models for alisim](http://www.iqtree.org/doc/Substitution-Models)
- -n - specifies the length of the root sequence (default: 500)
- -k - specifies how many simulations to be ran (default: 1)
- -s (optional) - can be true or false, if true it will also compute the saturation tests (default: false)

#### `analysis_biological.sh`
The analysis of real data requires the multiple sequence alignment file and optionally a tree file. If there is no tree file provided to program will ask to run IQ-Tree for the ML Tree and use that for further computations.

```sh
bash analysis_biological.sh -a alignment.phy -t treefile.tree -s true
```
- -a - specifies the sequence alignment file in PHYLIP or NEX format
- -t (optional) - the tree file in standard Newick format
- -s (optional) - can be true or false, if true it will also compute the saturation tests (default: false)

## Output
Both scripts produce up to 4 `.csv` files.
- `results_raw_<treename>.csv` - contains the raw values of the test statistics for each pair
- `results_raw_<treename>.csv` - contains the p-values against the null hypotheses with significance 0.05
- `results_rev_test.csv` - contains the results of the decision for each pair and test (whether the null hypothesis is retained/rejected)
- `results_sat_test.csv` - contains the results of the decision for each pair and each saturation test (if chosen)

Additionally one PDF file for each of the tests computed:
- `plot_Bowker_test.pdf` - coloured tree plot, heatmap and distribution of test statistics for Bowker test
- `plot_Stuart_test.pdf` - coloured tree plot, heatmap and distribution of test statistics for Stuart test
- `plot_IS_test.pdf` - coloured tree plot, heatmap and distribution of test statistics for Test for Internal Symmetry
- `plot_QS_test.pdf` - coloured tree plot, heatmap and distribution of test statistics for Test for Quasi-Symmetry
- `plot_Sat_Cassius1_test.pdf` - coloured tree plot, heatmap and distribution of test statistics for Cassius' Test for Saturation 1
- `plot_Sat_Cassius2_test.pdf` - coloured tree plot, heatmap and distribution of test statistics for Cassius' Test for Saturation 1
- `plot_Sat_Chi_test.pdf` - coloured tree plot, heatmap and distribution of test statistics for Chi-square Test for Saturation
- 'venn_diag.pdf' - Venn diagram for each pair, comparing the number of rejections in Bowker, Stuart and Quasi-Symmetry Test

And lastly the result tree with added labels on the branches for all test rejections in NEXUS format.
- `<treename>.tree`

## Example Input

To test if everything is working correctly, there is a test input provided with an example tree file and example alignment. To run the software open a terminal and navigate to the directory, then execute the following commands:

### Simulation Study

For the simulation study:
```
bash analysis_simulation.sh -t test_input/example-treefile.tree -m JC -n 1000 -k 100 -s true
```
If everything worked there should be a new folder called `results_example-treefile`, containing all of the outputs.

### Biological Study

For the biological study:
```
bash analysis_biological.sh -a test_input/example-alignment.phy -t example-treefile.tree -s true
```
If everything worked there should be a new folder called `results_example-alignment`, containing all of the outputs.
