__version__ = '0.0.1'

import random, pickle, time, warnings, copy, json
import numpy as np
import pandas as pd
from Bio import SeqIO, SeqRecord, Align
from sequence_preprocessing import concatenate_fastq, read_fastq
from itertools import compress

# define an "aligner" to run DNA alignment algorithm
aligner = Align.PairwiseAligner()
aligner.mode = 'global' # cf. Needleman-Wunsch algorithm
aligner.match_score = 1
aligner.mismatch_score = -1
aligner.gap_score = -2

# specific scores that I use for primer alignment (consensus team only)
aligner.query_end_gap_score = 0 # must be 0 (to allow gap around the primer for global alignment) [query ~ sequence here]
aligner.target_end_gap_score = -100 # set to -100 if the primer should be fully integrated in the sequence [target ~ primer here]

def check_primer_in_seq(sequences, primer):
    """
    returns an array (of the same shape of sequences) of booleans indicating if the primer is contained in the sequences
    """
    norm_align_scores = np.zeros(len(sequences))
    for sequence_id in range(len(sequences)):
        sequence = sequences[sequence_id]

        # align a primer on a sequence
        alignments = aligner.align(sequence, primer)
        best_alignment = max(alignments, key=lambda x: x.score)
        alignment_score = best_alignment.score
        alignment_score_normalized = alignment_score/len(primer)
        norm_align_scores[sequence_id] = alignment_score_normalized

    threshold = 0.45 #TODO figure out how to set up the threshold, do some statistical tests?
    norm_align_scores[norm_align_scores > threshold] = 1
    norm_align_scores[norm_align_scores <= threshold] = 0 
    return norm_align_scores.astype(int)

def analyze_fastq_statistics(fastq_filepath, dna_reference_sequence):
    dna_reference_sequence = ''.join([char for char in dna_reference_sequence if char in "ATGC"]) # remove space, line break...
    sequences = read_fastq(fastq_filepath)
    sequences_str = [str(s.seq) for s in sequences]
    analysis = dict()

    # analysis on sequences
    analysis['general'] = dict()
    analysis['general']['n_sequences'] = len(sequences_str)
    analysis['general']['average_sequence_length'] = np.mean([len(s) for s in sequences_str])
    
    forward_primer = 'ATGTCACCACAAACAGAGACTAAAGC'
    reverse_primer = 'TCGCATGTACCTGCAGTAGC'
    
    has_forward_primer = check_primer_in_seq(sequences_str, forward_primer)
    has_reverse_primer = check_primer_in_seq(sequences_str, reverse_primer)
    num_with_forward_primer = np.sum(has_forward_primer)
    ratio_forward_primer = np.mean(has_forward_primer)
    num_with_reverse_primer = np.sum(has_reverse_primer)
    ratio_reverse_primer = np.mean(has_reverse_primer)
    has_both_primers = np.logical_and(has_forward_primer, has_reverse_primer) - 0

    has_no_primer = 1 - has_both_primers
    num_with_both_primers = np.count_nonzero(has_both_primers)
    ratio_both_primers = np.mean(has_both_primers)
    num_with_no_primer = np.count_nonzero(has_no_primer)
    ratio_no_primer = np.mean(has_no_primer)

    analysis['general']['n_sequences_with_f_primer'] = num_with_forward_primer
    analysis['general']['n_sequences_with_r_primer'] = num_with_reverse_primer
    analysis['general']['n_sequences_with_boths_primers'] = num_with_both_primers
    analysis['general']['n_sequences_with_no_primer'] = num_with_no_primer

    seqs_with_boths_primers = np.array(sequences_str)[has_both_primers] #TODO why with list and not np?

    return analysis


if __name__ == '__main__':
    tomato_rbcL_reference_sequence = (
        """
        ATGTCACCACAAACAGAGACTAAAGCAAGTGTTGGATTCAAAGCTGGTGT
        TAAAGAGTACAAATTGACTTATTATACTCCTGAGTACCAAACCAAGGATA
        CTGATATATTAGCAGCATTCCGAGTAACTCCTCAACCTGGAGTTCCACCT
        GAAGAAGCAGGGGCCGCGGTAGCTGCCGAATCTTCTACTGGTACATGGAC
        AACTGTATGGACCGATGGACTTACCAGTCTTGATCGTTACAAAGGGCGAT
        GCTACCGCATCGAGCGCGTTGTTGGAGAAAAAGATCAATATATTGCTTAT
        GTAGCTTACCCTTTAGACCTTTTTGAAGAAGGTTCCGTTACCAATATGTT
        TACTTCCATTGTAGGTAACGTATTTGGGTTCAAAGCCCTGCGCGCTCTAC
        GTCTGGAAGATCTGCGAATCCCTCCTGCTTATGTTAAAACTTTCCAAGGT
        CCGCCTCATGGGATCCAAGTTGAAAGAGATAAATTGAACAAGTATGGTCG
        TCCCCTGTTGGGATGTACTATTAAACCTAAATTGGGGTTATCTGCAAAAA
        ACTACGGTAGAGCTGTTTATGAATGTCTTCGCGGTGGACTTGATTTTACC
        AAAGATGATGAGAACGTGAACTCACAACCATTTATGCGTTGGAGAGATCG
        TTTCTTATTTTGTGCCGAAGCACTTTTTAAAGCACAGACTGAAACAGGTG
        AAATCAAAGGGCATTACTTGAATGCTACTGCAGGTACATGCGA
        """
    )

    analysis = analyze_fastq_statistics(fastq_filepath="rbcL_Qiagen_tomato_5000.fastq", dna_reference_sequence=tomato_rbcL_reference_sequence)
    print(json.dumps(analysis, indent=4)) # "pretty print" to visualize a dict TODO understand error

