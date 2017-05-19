# Computational Biology Code

This repository contains all the work and code I did for the University of
Southampton's Computational Biology module.

## proteinCompare.py
[`proteinCompare.py`](https://github.com/aloisklink/COMP3212CompBiologyLabs/blob/master/proteinCompare.py)
contains functions to:
* load data from the [FASTA format](https://en.wikipedia.org/wiki/FASTA_format),
* aligning two sequences using the [Smith–Waterman algoirithm](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm) or [Needleman-Wunsch algoirithm](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm).
* using a [hidden Markov model (HMM)](https://en.wikipedia.org/wiki/Hidden_Markov_model) and the [Viterbi algorithm](https://en.wikipedia.org/wiki/Viterbi_algorithm) to find the most likely sequence of states in a given sequence

### Requirements
Python 3 (should work in Python 2 as well but don't count on it)

## simulateChemicals.py
[`simulateChemicals.py`](https://github.com/aloisklink/COMP3212CompBiologyLabs/blob/master/simulateChemicals.py) 
contains functions to simulate a system of chemical equations 
(e.g. ∅ → X at a rate 1, X → Y at rate 2).

It can simulate the system both by using the ODE format, 
or by using [Gillespie's algorithm](https://en.wikipedia.org/wiki/Gillespie_algorithm).

### Requirements
Python 3, numpy, scipy, and pylab.
These all come with [Anaconda](https://www.continuum.io/downloads).

## secondaryStruct.py
[`secondaryStruct.py`](https://github.com/aloisklink/COMP3212CompBiologyLabs/blob/master/secondaryStruct.py)
contains functions for predicting the [secondary protein structures](https://en.wikipedia.org/wiki/Protein_secondary_structure) of a protein using a [multi-layer perceptron neural network](https://en.wikipedia.org/wiki/Multilayer_perceptron).
In addition, it also has a few utility functions.

### Requirements
Python 3, numpy, and [Tensorflow](https://www.tensorflow.org/install/)

# Labs
Code for COMP3212 Computational Biology Labs.
See the [COMP3212 website](https://secure.ecs.soton.ac.uk/notes/comp3212/tutorial/) for more details.

## Lab 1
Solution is in [`lab1.py`](https://github.com/aloisklink/COMP3212CompBiologyLabs/blob/master/lab1.py).

1. Write a program to implement Needleman-Wunsch for proteins
  * You will need the [blosum50](https://secure.ecs.soton.ac.uk/notes/comp3212/tutorial/blosum50.txt) scoring matrix
  * You can use any programming language
  * Run this on HEAGAWGHEE versus PAWHEAE
  * Compare this to page 23 in [lecture 5][lecture_5]
  * Match the protein sequence PQPTTPVSSFTSGSMLGRTDTALTNTYSAL with PSPTMEAVEASTASHPHSTSSYFATTYYHL 
2. Modify your program to implement the Smith-Waterman algorithm
  * Again run this on HEAGAWGHEE versus PAWHEAE
  * Compare this to page 5 in [lecture 6][lecture_6]
  * Find the best local match between MQNSHSGVNQLGGVFVNGRPLPDSTRQKIVELAHSGARPCDISRILQVSNGCVSKILGRY and TDDECHSGVNQLGGVFVGGRPLPDSTRQKIVELAHSGARPCDISRI

[lecture_5]: https://secure.ecs.soton.ac.uk/notes/comp3212/lectures/apb-lectures/hmm.pdf
[lecture_6]: https://secure.ecs.soton.ac.uk/notes/comp3212/lectures/apb-lectures/localMatching.pdf

3. We are going to test the BLAST algorithm
  * Download the Pax6 protein for the moust by going to <http://www.uniprot.org/uniprot/P63015> choose the "Format" tab and choose the FASTA (canonical) format
  * Do the sacme for the eyeless protein for the fruit fly <http://www.uniprot.org/uniprot/O96791>
  * Perform a BLAST sequence comparison using the web service at <https://blast.ncbi.nlm.nih.gov>
  
4. Program the following HMM for detecting CG rich regions  
![Hidden Markov Model]
(http://users.ecs.soton.ac.uk/ak9g14/COMP3212/hmmCG.png)

5. Write a viterbi algorithm for finding the most likely CG regions and find a way of drawing this

6. Run this on the genome for the [phase lambda](https://www.ncbi.nlm.nih.gov/nuccore/215104?report=fasta) or [here](https://secure.ecs.soton.ac.uk/notes/comp3212/tutorial/phaseLambda.fasta)

## Lab 2
Solution is in [`lab2/`](https://github.com/aloisklink/COMP3212CompBiologyLabs/tree/master/lab2)[`lab2.m`](https://github.com/aloisklink/COMP3212CompBiologyLabs/blob/master/lab2/lab2.m).

Requires MATLAB R2016b or later, as there is a [function inside the script](https://www.mathworks.com/help/matlab/matlab_prog/local-functions-in-scripts.html).

1. The Fisher-Wright Model
  * This tutorial asks you to explore the Fisher-Wright model where we consider a single locus with two alleles. The population size remains fixed. The mutant has a selection pressure of (1+s) meaning that the mutant produces (1+s) offspring for each offspring of the non-mutant. The forward mutation rate measures the chance of an non-mutate mutating to a mutant. The backward mutation rate is the probability of the reverse process.
  * Experiment with the `simulation()` function and find out when doesn't the mutant take over the population 
2. Markov Chain Analysis
  * Show that the matrix created by `transition_matrix()` is stochastic
  * Generate the probabilities distributions p(t)
  * Compute the steady state distribution for `P=100, s=0.1, v=0.01, u=0.01`
3. Diffusion Analysis
  * Experiment with the quality of the diffusion approximation
  * Find out what the diffusion equation gives you which you can't get from the Markov analysis (Markov Analysis becomes super slow as P increases)
4. Real Populations
  * Compare the functions `ga()` and `multisim()`

## Lab 3
Solution is in [`lab3.py`](https://github.com/aloisklink/COMP3212CompBiologyLabs/blob/master/lab3.py).

[Pweave](http://mpastell.com/pweave/) can be used to publish this Python script
with the following command:
 * `pypublish FIR_design.py` to make a html file
 * `pypublish -f pdf FIR_design.py` to make a pdf file

1. Write down the differential equation describing the system of chemical equations (assuming a volume of 1)
  * ∅ → X at a rate 1
  * X → Y at rate 2
  * 2 X + Y → 3 X at rate 0.02
  * X → ∅ at rate 0.04
2. Use a package to solve the differential equation for 500 time units starting from X(0)=Y(0)=0 (matlab will do this)
3. Write a Gillespie algorithm to simulate the same four chemical equation and plot the results for 500 time units (note that this is a lot of data to plot and you might want to save and plot the data only after X or Y have changed in number by at least 5.

## Major Assignment
The report can be found in [`secondaryStructAnalysis.Pnw`](https://github.com/aloisklink/COMP3212CompBiologyLabs/blob/master/secondaryStructAnalysis.Pnw).

It can be run with the following commands:

```bash
pweave -f texminted secondaryStructAnalysis.Pnw
pdflatex -synctex=1 -shell-escape secondaryStructAnalysis.tex
bibtex secondaryStructAnalysis.aux
pdflatex -synctex=1 -shell-escape secondaryStructAnalysis.tex
pdflatex -synctex=1 -shell-escape secondaryStructAnalysis.tex
```

The assignment is about secondary protein structure prediction. There are two main tasks:

1. Find, at least one, secondary protein structure prediction algorithm on the web and test is performance
2. Build you own secondary protein structure prediction algorithm and test its performance

The administrative details of the assessment are as follows
* Write a report of no more than 5 pages in at least 11pt on A4
* A PDF version of your report be submitted submitted to C-BASS by 4:00pm on Thursday 18th May 2017 (the penultimate day of term)
* This is worth 40% of the marks and should take you around 60 hours
* You are expected to work individually

In marking the following will be taken into consideration:
* You can use any test and training set you wish (for example UCI). You should briefly describe the data set you have chosen
* For the algorithm(s) you take from the web you should read the associated paper and give a brief summary of how it works
* The performance of your algorithm is not important
* You can use any machine learning algorithm you like (e.g. MLP, SVM, HMM)
* You can work in any programming language you want
* You may replicate an approach taken in a paper provided you re-implement it or you can follow your own approach
* You can use any standard package to perform machine learning (or write you own)
* You should acknowledge any packages you use
* You are expected to tackle the problem of representing the data yourself (that is, your package should not be tailored to secondary protein structure prediction)
* Marks will be awarded for testing different approaches, reflecting on the algorithm you have created, and comparing it to the algorithm(s) you have downloaded

