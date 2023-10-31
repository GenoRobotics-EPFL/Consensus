__version__ = '0.0.1'

import random, pickle, time, warnings, copy, json
import numpy as np
import pandas as pd
from Bio import SeqIO, SeqRecord, Align
from sequence_preprocessing import concatenate_fastq, read_fastq

def check_primer_in_seq(sequences, primer):
    """
    returns an array (of the same shape of sequences) of booleans indicating if the primer is contained in the sequences
    """
    has_primer = np.zeros(len(sequences))
    for s_idx, s in enumerate(sequences):
        if primer in s:
            has_primer[s_idx] = 1
    return has_primer

def analyze_fastq_statistics(fastq_filepath, dna_reference_sequence, dna_reference_genome_description):
    dna_reference_sequence = ''.join([char for char in dna_reference_sequence if char in "ATGC"]) # remove space, line break...
    sequences = read_fastq(fastq_filepath)
    sequences_str = [str(s.seq) for s in sequences]
    analysis = dict()

    # analysis on sequences
    analysis['general'] = dict()
    analysis['general']['n_sequences'] = len(sequences_str)
    analysis['general']['average_sequence_length'] = np.mean([len(s) for s in sequences_str])
    
    # TODO: figure out how to get the primers, is there already a function for it?
    forward_primer = ...
    reverse_primer = ...
    
    has_forward_primer = check_primer_in_seq(sequences_str, forward_primer)
    has_reverse_primer = check_primer_in_seq(sequences_str, reverse_primer)
    num_with_forward_primer = np.sum(has_forward_primer)
    ratio_forward_primer = np.mean(has_forward_primer)
    num_with_reverse_primer = np.sum(has_reverse_primer)
    ratio_reverse_primer = np.mean(has_reverse_primer)
    has_both_primers = np.logical_and(has_forward_primer, has_reverse_primer)
    has_no_primer = 1 - has_both_primers
    num_with_both_primers = np.sum(has_both_primers)
    ratio_both_primers = np.mean(has_both_primers)
    num_with_no_primer = np.sum(has_no_primer)
    ratio_no_primer = np.mean(has_no_primer)

    analysis['general']['n_sequences_with_f_primer'] = num_with_forward_primer
    analysis['general']['n_sequences_with_r_primer'] = num_with_reverse_primer
    analysis['general']['n_sequences_with_boths_primers'] = num_with_both_primers
    analysis['general']['n_sequences_with_no_primer'] = num_with_no_primer

    seqs_with_boths_primers = sequences[has_both_primers]

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
    print(json.dumps(analysis, indent=4)) # "pretty print" to visualize a dict

