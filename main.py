# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch
import numpy as np

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    
    my_NW = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    hs_gg = my_NW.align(hs_seq, gg_seq)[0]
    hs_mm = my_NW.align(hs_seq, mm_seq)[0]
    hs_br = my_NW.align(hs_seq, br_seq)[0]
    hs_tt = my_NW.align(hs_seq, tt_seq)[0]

    sim_scores = [(hs_gg, 'Gallus_gallus'),
                   (hs_mm, 'Mus_musculus'),
                   (hs_br, 'Balaeniceps_rex'),
                   (hs_tt, 'Tursiops_truncatus')]
    sorted_by_sim = sorted(sim_scores, key=lambda tup: tup[0], reverse=True)

    print("Species in order of most similar to least similar: ")
    print([x[1] for x in sorted_by_sim])

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    
    print("Alignment scores between each species BRD2 and human BRD2 (score, non-human species): ")
    print(sorted_by_sim)


if __name__ == "__main__":
    main()
