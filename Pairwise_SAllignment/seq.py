#Match/Gap Parameters Descriptions
#x     No parameters. Identical characters have score of 1, otherwise 0.
#m     A match score is the score of identical chars, otherwise mismatch
#      score.
#d     A dictionary returns the score of any pair of characters.
#c     A callback function returns scores.

#x     No gap penalties.
#s     Same open and extend gap penalties for both sequences.
#d     The sequences have different open and extend gap penalties.
#c     A callback function returns the gap penalties.

##Global Allignment

#from Bio import pairwise2
#alignments = pairwise2.align.globalxx("ACCGT", "ACG")
#print(alignments)

#Local Allignment
from Bio import SubsMat
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
alignments = pairwise2.align.globalms("ATG", "AGG", 2, -1, -1, 0)  # Uses the 'm' scoring scheme
print(format_alignment(*alignments[0]))

for a in pairwise2.align.globalxx("ACCGT", "ACG"):

    print(format_alignment(*a))

for a in pairwise2.align.localxx("ACCGT", "ACG"):

    print(format_alignment(*a))

for a in pairwise2.align.localxx("ACCGT", "ACG"):

    print(format_alignment(*a, full_sequences=True))

from Bio.SubsMat import MatrixInfo as matlist

matrix = matlist.blosum62

for a in pairwise2.align.globaldx("KEVLA", "EVL", matrix):

    print(format_alignment(*a))

from math import log

def gap_function(x, y):  # x is gap position in seq, y is gap length

    if y == 0:  # No gap

        return 0

    elif y == 1:  # Gap open penalty

        return -2

    return - (2 + y/4.0 + log(y)/2.0)


alignment = pairwise2.align.globalmc("ACCCCCGT", "ACG", 5, -4,

                                     gap_function, gap_function)