

import gzip, glob, math, random, pickle, time, warnings, sys, copy
import numpy as np
from Bio import SeqIO, SeqRecord, Align
import sequence_preprocessing as sp

class PrimerAlignment():
   
    def __init__(self, primer_score=None, primer_start=None, primer_end=None, primer_confidence=None, primer_alignment=None, sequence_complementary=False, sequence_reverse=False):
        self.primer_score = primer_score  # alignment score between primer and sequence [float]
        self.primer_start = primer_start  # start position of primer in sequence [int]
        self.primer_end = primer_end  # end position of primer in sequence [int]
        self.primer_confidence = primer_confidence
        self.primer_alignment = primer_alignment  # alignment information formatted (between primer and sequence) [list of str]
        self.sequence_complementary = sequence_complementary
        self.sequence_reverse = sequence_reverse
   
   
    
    def align_primer(sequence, primer, aligner_scores=None, print_alignment=False, return_formatted_alignment=False):
        """
        Global alignment of primer inside sequence.
        :param sequence: sequence where the primer is searched [str]
        :param primer: primer searched [str]
        :param print_alignment: print alignment results [bool]
        :return: primer alignment results [float, int, int, [str, str, str]]
                score: alignment score [float]
                loc_start: start position of the primer in the sequence [int]
                loc_end: end position of the primer in the sequence [int]
                primer_alignment: list of 3 strings containing formatted sequence alignment (truncated sequence, alignment, primer) [list of 3 str]
        """
        # Aligner settings
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        # print(Align.substitution_matrices.load()) # print possible substitution matrices
        if aligner_scores is None: aligner_scores = dict()
        aligner.match_score = aligner_scores.get("match_score") if aligner_scores.get("match_score") is not None else 1 #3
        aligner.mismatch_score = aligner_scores.get("mismatch_score") if aligner_scores.get("mismatch_score") is not None else -1 #-2
        aligner.gap_score = aligner_scores.get("gap_score") if aligner_scores.get("gap_score") is not None else -2
        aligner.query_end_gap_score = 0 # must be 0 (to allow gap around the primer for global alignment) [query ~ sequence here]
        aligner.target_end_gap_score = -100 # set to -100 if the primer should be fully integrated in the sequence [target ~ primer here]

        # Alignment
        alignments = aligner.align(sequence, primer)
        # Get the best alignment only
        alignment = max(alignments, key=lambda x: x.score) # old version (slower): # alignment = sorted(alignments)[0]
        if print_alignment:
            print(f"Score = {round(alignment.score,1)}")
            print(alignment)

        # 'score', 'loc_start', 'loc_end' computations
        score = alignment.score/len(primer)
        loc_start = alignment.aligned[0][0][0]
        loc_end = alignment.aligned[0][-1][1]

        # 'primer_alignment' contains formatted string used to print the alignment results
        primer_alignment = None
        if return_formatted_alignment:
            primer_alignment = []
            aform = alignment.format() # format() calls _format_pretty()
            aform_linebreak1 = aform.find('\n')
            aform_linebreak2 = aform.find('\n', aform_linebreak1 + 1)
            aform_linebreak3 = aform.find('\n', aform_linebreak2 + 1)
            primer_position_str = aform[aform_linebreak2+1:aform_linebreak3].replace("A", "N").replace("T", "N").replace("G", "N").replace("C", "N")
            pos1 = primer_position_str.find("N")
            pos2 = primer_position_str.rfind("N") + 1
            primer_alignment.append(aform[0:aform_linebreak1][pos1:pos2])
            primer_alignment.append(aform[aform_linebreak1+1:aform_linebreak2][pos1:pos2])
            primer_alignment.append(aform[aform_linebreak2+1:aform_linebreak3][pos1:pos2])

        conf = sp.sigmoid_primer(score)

        return PrimerAlignment(primer_score=score, primer_start=loc_start, primer_end=loc_end, primer_confidence=conf, primer_alignment=primer_alignment)
    
    def get_formatted_alignment(self, sequence, primer):
        """
        Get formatted strings of the primer alignment
        :param sequence: sequence [str]
                         the sequence given must be "reverse" if self.sequence_reverse = True
                         the sequence given must be "complementary" if self.sequence_complementary = True
        :return: list of 3 strings containing formatted sequence alignment (truncated sequence, alignment, primer) [list of 3 str]
        """
        if self.sequence_reverse: sequence = sequence[::-1]
        if self.sequence_complementary: sequence = sequence.replace("A", "t").replace("T", "a").replace("G", "c").replace("C", "g").upper()
        seq_primed = sequence[self.primer_start:self.primer_end]

        return PrimerAlignment.align_primer(seq_primed, primer, return_formatted_alignment=True).primer_alignment