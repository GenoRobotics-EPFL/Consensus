__version__ = '0.0.1'

import random, pickle, time, warnings, copy, json
import numpy as np
import pandas as pd
from Bio import SeqIO, SeqRecord, Align
from sequence_preprocessing import concatenate_fastq, read_fastq


def analyze_fastq_statistics(fastq_filepath, dna_reference_sequence):
    dna_reference_sequence = ''.join([char for char in dna_reference_sequence if char in "ATGC"]) # remove space, line break...
    sequences = read_fastq(fastq_filepath)
    sequences_str = [str(s.seq) for s in sequences]
    analysis = dict()

    # analysis on sequences
    analysis['general'] = dict()
    analysis['general']['n_sequences'] = len(sequences_str)
    analysis['general']['average_sequence_length'] = np.mean([len(s) for s in sequences_str])

    # TODO: add other analysis below :)


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

