# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")

    M_truth = np.array([[0, -np.inf, -np.inf, -np.inf, -np.inf], 
                      [-np.inf, 5, -12, -12, -14],
                      [-np.inf, -11, 4, -1, -6],
                      [-np.inf, -13, -8, 5, 4]])
    X_truth = np.array([[-10, -np.inf, -np.inf, -np.inf, -np.inf], 
                      [-11, -22, -23, -24, -25],
                      [-12, -6, -17, -18, -19],
                      [-13, -7, -7, -12, -17]])
    Y_truth = np.array([[-10, -11, -12, -13, -14], 
                      [-np.inf, -22, -6, -7, -8],
                      [-np.inf, -23, -17, -7, -8],
                      [-np.inf, -24, -18, -18, -6]])

    my_NW = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    my_NW.align(seq1, seq2)

    assert np.array_equal(my_NW._align_matrix, M_truth)
    assert np.array_equal(my_NW._gapA_matrix, X_truth)
    assert np.array_equal(my_NW._gapB_matrix, Y_truth)

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    
    my_NW = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    backtrace = my_NW.align(seq3, seq4)

    assert backtrace == (17, "MAVHQLIRRP", "M---QLIRHP")





