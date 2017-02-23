# COMP3212CompBiologyLabs
Code for COMP3212 Computational Biology Labs.
See the [COMP3212 website](https://secure.ecs.soton.ac.uk/notes/comp3212/tutorial/) for more details.

Mirror below for convience:

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
