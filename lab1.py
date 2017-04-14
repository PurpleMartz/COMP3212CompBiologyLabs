from proteinCompare import *

blosumCosts = loadBlosum50()
match = AminoAcidMutation(blosumCosts)
print("Needleman-Wunsch: ")
match.needlemanWunsch("HEAGAWGHEE", "PAWHEAE")
match.needlemanWunsch("PQPTTPVSSFTSGSMLGRTDTALTNTYSAL",
                      "PSPTMEAVEASTASHPHSTSSYFATTYYHL")
print("")
print("Smith-Waterman: ")
match.smithWaterman("HEAGAWGHEE", "PAWHEAE")
match.smithWaterman("MQNSHSGVNQLGGVFVNGRPLPDSTRQKIVELAHSGARPCDISRILQVSNGCVSKILGRY",
                    "TDDECHSGVNQLGGVFVGGRPLPDSTRQKIVELAHSGARPCDISRI")
print("")
print("Hidden Markov Model")
dishonestCasino = HiddenMarkovModel(("123456"), 
                        (HiddenMarkovModel.State((1/6, 1/6, 1/6, 1/6, 1/6, 1/6),
                                                 (9/10, 1/10)),
                         HiddenMarkovModel.State((1/10, 1/10, 1/10, 1/10, 1/10, 1/2),
                                                 (1/10, 9/10))))
dishonestCasino.viterbi("5453525456666664365666635661416626365666621166211311155566351166565663466653642535666662541345464155")

cgRich = HiddenMarkovModel(("ATCG"),
                     (HiddenMarkovModel.State((0.2698, 0.3237, 0.2080, 0.1985),
                                              (0.9998, 0.0002)),
                     (HiddenMarkovModel.State((0.2459, 0.2079, 0.2478, 0.2984),
                                              (0.0003, 0.9997)))))
cgRich.viterbi(Fasta("./data/phaseLambda.fasta").sequence)
