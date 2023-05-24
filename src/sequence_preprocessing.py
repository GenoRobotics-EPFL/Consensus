"""
This version of the file is a TEST version.
It means that some methods are not finished or not optimized.
It also means that everything here is subject to change in later versions.
- Jeremy Goumaz
"""

__version__ = '2.1.0'

import gzip, glob, math, random, pickle, time, warnings, sys, copy
import numpy as np
from Bio import SeqIO, SeqRecord, Align



def read_fastq(fastq_filepath=None):
    """
    Read a fastq file and return the reads.
    :param fastq_filepath: filepath of the .fastq file [str]
    :return: reads from the fastq file [list of Bio.SeqRecord.SeqRecord]
    """
    if fastq_filepath is None: fastq_filepath = "data/rbcL_Qiagen_tomato.fastq" # default path (example)
    if fastq_filepath.lower().endswith('.gz'):
        f = gzip.open(fastq_filepath, 'rt')
    else:
        f = open(fastq_filepath, 'rt')
    reads = []
    for read in SeqIO.parse(f, "fastq"):
        reads.append(read)
    return reads
def concatenate_fastq(src_folder=None, dst=None):
    """
    Concatenate all .fastq from a folder into a single .fastq file.
    :param folder: folder containing the .fastq files (usually fastq_pass folder) [str]
    :param dst: destination file (.fastq) [str]
    :return: None
    """
    if src_folder is None: src_folder = "fastq_pass" # default folder (example)
    if dst is None: dst = f"{src_folder}/concatenation.fastq" # default destination (example)
    if dst[-6:] != ".fastq": dst = f"{dst}.fastq"
    def get_file_iterator(filename):
        if filename.lower().endswith('.gz'):
            f = gzip.open(filename, 'rt')
        else:
            f = open(filename, 'rt')
        return f
    fastq_list = [fastq for fastq in glob.glob(f"{src_folder}/*.fastq*")]
    # print(fastq_list)
    fastq_iterators = [SeqIO.parse(get_file_iterator(fastq), "fastq") for fastq in fastq_list]
    while True:
        for fq in fastq_iterators:
            try:
                # print(next(fq).format("fastq"), end="")
                SeqIO.write(next(fq), open(dst,"at"), "fastq")
            except StopIteration:
                fastq_iterators.remove(fq)
        if len(fastq_iterators) == 0:
            break
def get_substring(string, k): # ACTUALLY NOT USED IN THE CODE
    """
    Get substring of length k from a string.
    :param string: original string [str]
    :param k: length of substrings [int]
    :return: list of substrings [list of str | len=len(string)-k+1]
    """
    return [string[i:i+k] for i in range(0, len(string)-k+1)]
def get_subprimer(primer, k, error_ratio=0): # ACTUALLY NOT USED IN THE CODE
    """
    Get subprimers (substrings of a primer) of length k and the
    :param primer: primer sequence [str]
    :param k: length of subprimers [int]
    :param error_ratio: error_ratio allowed [float]
    :return: list of subprimers (subprimer sequence, left characters to try, right characters to try) [list of [str, int, int]]
    """
    ratio = 1.0 + error_ratio
    return [[primer[i:i+k], math.ceil(ratio*i), math.ceil(ratio*(len(primer)-i-k))] for i in range(0, len(primer)-k+1)]
def create_random_sequence(length=500, seed=None):
    """
    Generate a random sequence.
    :param length: length of the sequence generated [int]
    :param seed: random seed [int]
    :return: random sequence [str]
    """
    if seed is not None: random.seed(seed)
    return "".join(["ATGC"[random.randint(0,3)] for i in range(length)])
def damage_sequence(sequence, mutation_rate=0.05, deletion_rate=0.05, insertion_rate=0.05):
    """
    Damage a sequence randomly with mutation (substitution), deletion and insertion
    Warning: to avoid infinite loop, don't set the insertion_rate to 1.0 (or near)
    :param sequence: original sequence to damage [str]
    :param mutation_rate: mutation rate (between 0.0 and 1.0) [float]
    :param deletion_rate: deletion_rate (between 0.0 and 1.0) [float]
    :param insertion_rate: insertion_rate (between 0.0 and 1.0) [float]
    :return: damaged sequence [str]
    """
    if mutation_rate < 0 or mutation_rate > 1:
        raise Exception("[damage_sequence] mutation_rate is incorrect (must be between 0.0 and 1.0)")
    if deletion_rate < 0 or deletion_rate > 1:
        raise Exception("[damage_sequence] deletion_rate is incorrect (must be between 0.0 and 1.0)")
    if insertion_rate < 0 or insertion_rate > 1:
        raise Exception("[damage_sequence] insertion_rate is incorrect (must be between 0.0 and 1.0)")
    sequence = "".join(["ATGC".replace(b, '')[random.randint(0, 2)] if random.random() < mutation_rate else b for b in sequence]) # mutation / substitution
    sequence = "".join(['' if random.random() < deletion_rate else b for b in sequence]) # deletion
    insertion_extension_rate = insertion_rate  # can be changed if the extension rate is different
    def get_insert(extension_rate=insertion_extension_rate):
        insert = "ATGC"[random.randint(0,3)]
        while random.random() < extension_rate: insert += "ATGC"[random.randint(0,3)] # extension (after insertion)
        return insert
    sequence = "".join([sequence[i:i+1] + get_insert() if random.random() < insertion_rate else sequence[i:i+1] for i in range(len(sequence)+1)]) # insertion
    return sequence
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

    conf = sigmoid_primer(score)
    return PrimerAlignment(primer_score=score, primer_start=loc_start, primer_end=loc_end, primer_confidence=conf, primer_alignment=primer_alignment)
def align_primer_random_scores(seq_length=500, primer_length=26, iterations=1000, seed=0, print_scores=True):
    """
    Evaluate alignment scores of random alignments (useful to determine score threshold to avoid the use of random alignments)
    :param seq_length: length of sequence [int]
    :param primer_length: length of primer [int]
    :param iterations: number of iterations (more iterations -> more statistically significant but slower) [int]
    :param seed: random seed [int]
    :param print_scores: print the results/scores [bool]
    :return: 99.9% percentile score [float]
    """
    if seed is not None: random.seed(seed)
    scores = np.empty(shape=(iterations,))
    for it in range(iterations):
        primer = create_random_sequence(primer_length)
        seq = create_random_sequence(seq_length)
        scores[it] = align_primer(seq, primer).primer_score
        # if it%1000==0: print(it)
    scores = scores/np.mean(scores)
    if print_scores:
        print(f"Random alignments scores (seq_length={seq_length}, primer_length={primer_length}, it={iterations}):")
        print(f"\tMean: {np.mean(scores)}")
        print(f"\tStd: {np.std(scores)}")
        for percent in [50, 75, 90, 95, 99, 99.9]:
            print(f"\tBest {round(100-percent,1)}%: {round(np.percentile(scores, percent),2)}")
    return scores
def sigmoid_primer(x):
    # temporary function (ASUP later) -> used to give a probability to a primer score
    x = 14 * x - 7
    if x >= 100:
        sig = 1.0
    elif x <= -100:
        sig = 0.0
    else:
        sig = 1/(1+np.exp(-x))
    return sig

class PrimerAlignment():
    def __init__(self, primer_score=None, primer_start=None, primer_end=None, primer_confidence=None, primer_alignment=None, sequence_complementary=False, sequence_reverse=False):
        self.primer_score = primer_score  # alignment score between primer and sequence [float]
        self.primer_start = primer_start  # start position of primer in sequence [int]
        self.primer_end = primer_end  # end position of primer in sequence [int]
        self.primer_confidence = primer_confidence
        self.primer_alignment = primer_alignment  # alignment information formatted (between primer and sequence) [list of str]
        self.sequence_complementary = sequence_complementary
        self.sequence_reverse = sequence_reverse
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
        return align_primer(seq_primed, primer, return_formatted_alignment=True).primer_alignment
class PrimerAllAlignments():
    def __init__(self, primer):
        self.primer = primer  # primer sequence [str]
        self.alignments = [[None, None], [None, None]] # 4 instances of PrimerAlignment [[normal, reverse], [complementary, reverse_complementary]]
    def set_alignment(self, PrimerAlign, complementary=False, reverse=False):
        """
        Setter of a PrimerAlignment
        :param PrimerAlign: PrimerAlignment to set [PrimerAlignment instance]
        :param complementary: is the alignment for complementary sequence [bool]
        :param reversed: is the alignment for reversed sequence [bool]
        :return: None
        """
        PrimerAlign.sequence_complementary = complementary
        PrimerAlign.sequence_reverse = reverse
        self.alignments[int(complementary)][int(reverse)] = PrimerAlign
    def get_formatted_main_alignment(self, sequence):
        if self.alignments[0][0] is not None:
            main_primer_alignment = self.alignments[0][0].get_formatted_alignment(sequence, self.primer)
        else:
            main_primer_alignment = None
            warnings.warn("PrimerAllAlignments.get_formatted_main_alignment(): the main primer_alignment has not been set. It can't be accessed.")
        return main_primer_alignment
    def get_main_alignment_properties(self):
        if self.alignments[0][0] is not None:
            properties = {"score": self.alignments[0][0].primer_score,
                          "start": self.alignments[0][0].primer_start,
                          "end": self.alignments[0][0].primer_end,
                          "confidence": self.alignments[0][0].primer_confidence}
            return properties
        else:
            warnings.warn("PrimerAllAlignments.get_main_alignment_properties(): the main primer_alignment has not been set. It can't be accessed.")
            return dict()
    def flip_reverse(self):
        new_alignments = [[self.alignments[0][1], self.alignments[0][0]], [self.alignments[1][1], self.alignments[1][0]]]
        self.alignments = new_alignments
        for id1 in range(2):
            for id2 in range(2):
                self.alignments[id1][id2].sequence_reverse = not self.alignments[id1][id2].sequence_reverse
    def flip_complementary(self):
        new_alignments = [[self.alignments[1][0], self.alignments[1][1]], [self.alignments[0][0], self.alignments[0][1]]]
        self.alignments = new_alignments
        for id1 in range(2):
            for id2 in range(2):
                self.alignments[id1][id2].sequence_complementary = not self.alignments[id1][id2].sequence_complementary
    def truncate(self, start, end):
        for id1 in range(2):
            for id2 in range(2):
                if self.alignments[id1][id2].sequence_reverse:
                    start_local, end_local = end, start
                else:
                    start_local, end_local = start, end
                if start_local > self.alignments[id1][id2].primer_start:
                    self.alignments[id1][id2] = None
                else:
                    self.alignments[id1][id2].primer_start = self.alignments[id1][id2].primer_start - start_local
                    self.alignments[id1][id2].primer_end = self.alignments[id1][id2].primer_end - start_local
    def extend(self, start_added, end_added):
        for id1 in range(2):
            for id2 in range(2):
                if self.alignments[id1][id2].sequence_reverse:
                    start_added_local, end_added_local = end_added, start_added
                else:
                    start_added_local, end_added_local = start_added, end_added
                self.alignments[id1][id2].primer_start = self.alignments[id1][id2].primer_start + start_added_local
                self.alignments[id1][id2].primer_end = self.alignments[id1][id2].primer_end + start_added_local
    def get_best_direction(self):
        best = {"confidence": None, "complementary": None, "reverse": None}
        for is_complementary in range(2):
            for is_reverse in range(2):
                primer_confidence = self.alignments[is_complementary][is_reverse].primer_confidence
                if primer_confidence is not None:
                    if best["confidence"] is None:
                        best["confidence"] = 0.0
                    if primer_confidence > best["confidence"]:
                        best["confidence"] = primer_confidence
                        best["complementary"] = bool(is_complementary)
                        best["reverse"] = bool(is_reverse)
        return best
class Sequence():
    def __init__(self, seq=None):
        self.sequence = None  # sequence [str]
        self.phred_scores = None  # Phred scores [np.ndarray]
        self.sequence_description = None  # sequence description from fastq files [str]
        self.mean_error_rate = None  # mean error rate computed from Phred scores [float]
        if seq is not None: self.set_sequence(seq)
        self.primer_forward_alignments = None # PrimerAllAlignments instance
        self.primer_reverse_alignments = None # PrimerAllAlignments instance
        self.primer_aligned = {"primer": None, "confidence": 0.0} # information about the primer aligned ("forward" or "reverse") and the primer alignment confidence associated
        self.primer_alignments_run = False # boolean specifying if run_primer_alignments() has been called
    def reset(self):
        """
        Reset the Sequence instance.
        :return: None
        """
        self.sequence = None
        self.phred_scores = None
        self.sequence_description = None
        self.mean_error_rate = None
        self.primer_forward_alignments = None
        self.primer_reverse_alignments = None
        self.primer_alignments_run = False
    def set_phred_scores(self, phred_scores):
        """
        Set Phred scores
        :param phred_scores: phred_scores [np.ndarray | 1D]
        :return: None
        """
        if type(self.phred_scores) != np.ndarray:
            phred_scores = np.array(phred_scores).reshape((-1,))
        self.phred_scores = phred_scores
        self.compute_mean_error_rate()
    def get_sequence_reversed(self):
        """
        Get the reverse sequence (5'-3' <-> 3'-5')
        /!\ Warning: reverse sequence has not the same meaning as the reverse primer
        :return: reversed sequence [str]
        """
        return self.sequence[::-1]
    def get_sequence_complementary(self):
        """
        Get the complementary sequence
        :return: complementary sequence [str]
        """
        return self.sequence.replace("A", "t").replace("T", "a").replace("G", "c").replace("C", "g").upper()
    def get_sequence_variants(self):
        """
        Get the 4 variants of the sequence:
            - sequence (code 0|00)
            - complementary sequence (code 1|10)
            - reversed sequence (code 2|01)
            - reversed complementary sequence (code 3|11)
        :return: list of sequences [list of str | len=4]
        """
        sequence_comp = self.get_sequence_complementary()
        return [self.sequence, sequence_comp, self.sequence[::-1], sequence_comp[::-1]]
    def compute_mean_error_rate(self):
        """
        Compute and set the mean error rate (using the phred scores)
        :return: None
        """
        if self.phred_scores is not None:
            error_proba = 10 ** (-self.phred_scores / 10)
            self.mean_error_rate = error_proba.mean()
        else:
            self.mean_error_rate = None
    def flip_reverse(self):
        """
        Flip the Sequence to the reverse.
        :return: None
        """
        self.sequence = self.get_sequence_reversed()
        if self.phred_scores is not None:
            self.phred_scores = np.flip(self.phred_scores)
        if self.primer_forward_alignments is not None: self.primer_forward_alignments.flip_reverse()
        if self.primer_reverse_alignments is not None: self.primer_reverse_alignments.flip_reverse()
    def flip_complementary(self):
        """
        Flip the Sequence to the complementary.
        :return: None
        """
        self.sequence = self.get_sequence_complementary()
        if self.primer_forward_alignments is not None: self.primer_forward_alignments.flip_complementary()
        if self.primer_reverse_alignments is not None: self.primer_reverse_alignments.flip_complementary()
    def truncate(self, start, end=None):
        """
        Truncate the sequence.
        :param start: start of the truncation
        :param end: end of the truncation
        :return: None
        """
        if end is None: end = len(self.sequence)
        if start > len(self.sequence) or end > len(self.sequence):
            warnings.warn("Sequence.truncate(): invalid start/end given as arguments.")
        self.sequence = self.sequence[start:end]
        if self.phred_scores is not None:
            self.phred_scores = self.phred_scores[start:end]
        self.compute_mean_error_rate()
        if self.primer_forward_alignments is not None: self.primer_forward_alignments.truncate(start, end)
        if self.primer_reverse_alignments is not None: self.primer_reverse_alignments.truncate(start, end)
    def extend(self, start_added=0, end_added=0):
        """
        Extend the sequence with unknown nucleotides (N).
        :param start_added: number of nucleotides to add in the beginning [int]
        :param end_added: number of nucleotides to add at the end [int]
        :return: None
        """
        self.sequence = start_added*"N" + self.sequence + end_added*"N"
        if self.phred_scores is not None:
            self.phred_scores = np.concatenate((np.zeros(shape=(start_added,)), self.phred_scores, np.zeros(shape=(end_added,))))
        self.compute_mean_error_rate()
        if self.primer_forward_alignments is not None: self.primer_forward_alignments.extend(start_added, end_added)
        if self.primer_reverse_alignments is not None: self.primer_reverse_alignments.extend(start_added, end_added)
    def set_sequence(self, seq):
        """
        Set the sequence
        :param seq: sequence [Bio.SeqRecord.SeqRecord OR str]
        :return: None
        """
        if type(seq) is str:
            seq = ''.join([s for s in str(seq) if s in "ATGC"])
            if self.sequence is not None:
                found = False
                solvable = True
                for id, var in enumerate(self.get_sequence_variants()):
                    if seq in var:
                        if found: solvable = False # not possible if seq present in multiple var
                        found = True
                        loc_start = var.find(seq)
                        loc_end = loc_start + len(seq)
                        if var[loc_start+1:].find(seq) != -1: solvable = False # not possible if seq present multiple times in var
                        flip_complementary = id%2==1
                        flip_reverse = id//2==1
                        transformations_to_do = [flip_complementary, flip_reverse, loc_start, loc_end]
                if found == False or solvable == False:
                    self.reset()
                    self.sequence = seq
                else:
                    if transformations_to_do[0] == True: self.flip_complementary()
                    if transformations_to_do[1] == True: self.flip_reverse()
                    self.truncate(start=transformations_to_do[2], end=transformations_to_do[3])
            else:
                self.sequence = seq
        elif type(seq) is SeqRecord.SeqRecord:
            self.sequence = str(seq.seq)
            self.phred_scores = np.array(seq.letter_annotations['phred_quality'])
            self.sequence_description = str(seq.description)
        else:
            warnings.warn("The sequence type given to Sequence() class is not recognized.")
        self.compute_mean_error_rate()
    def set_sequence_with_best_primer_alignment(self):
        # Get the best directions from each primer alignments
        best_confidence_forward = self.primer_forward_alignments.get_best_direction()
        best_confidence_reverse = self.primer_reverse_alignments.get_best_direction()

        # Check if the best directions have been found
        error_forward, error_reverse = False, False
        if best_confidence_forward["confidence"] is None:
            best_confidence_forward["confidence"] = 0.0
            error_forward = True
        if best_confidence_reverse["confidence"] is None:
            best_confidence_reverse["confidence"] = 0.0
            error_reverse = True
        if error_forward and error_reverse:
            warnings.warn("[set_sequence_with_best_primer_alignments]: all primer alignments are empty")
            return None
        else:
            if best_confidence_forward["confidence"] > best_confidence_reverse["confidence"]:
                best_confidence = best_confidence_forward
                self.primer_aligned = {"primer": "forward", "confidence": best_confidence["confidence"]}
            else:
                best_confidence = best_confidence_reverse
                self.primer_aligned = {"primer": "reverse", "confidence": best_confidence["confidence"]}

        # Set the new sequence direction
        if best_confidence["complementary"]:
            self.flip_complementary()
        if best_confidence["reverse"]:
            self.flip_reverse()
    def truncate_sequence_using_best_primer_alignment(self):
        self.set_sequence_with_best_primer_alignment()
        if self.primer_aligned["primer"] == "forward":
            start = self.primer_forward_alignments.get_main_alignment_properties().get("start")
            self.truncate(start=start)
        elif self.primer_aligned["primer"] == "reverse":
            end = self.primer_reverse_alignments.get_main_alignment_properties().get("end")
            self.truncate(start=0, end=len(self.sequence)-end)
    def run_primer_alignments(self, primer_forward, primer_reverse):
        """
        Run primer_align for each Sequence.get_variants() for each primer (forward and reverse)
        This method initializes/sets the primers properties.
        :param primer_forward: 5'->3' forward primer [str]
                               ex: "ATGTCACCACAAACAGAGACTAAAGC"
        :param primer_reverse: 5'->3' reverse primer [str]
                               ex: "TCGCATGTACCTGCAGTAGC"
        :return: None
        """
        # Complementarize/reverse the primer_reverse (in order to put it in the same direction as primer_forward)
        primer_reverse = primer_reverse.replace("A", "t").replace("T", "a").replace("G", "c").replace("C", "g").upper()
        primer_reverse = primer_reverse[::-1]

        self.primer_forward_alignments = PrimerAllAlignments(primer_forward)
        self.primer_reverse_alignments = PrimerAllAlignments(primer_reverse)
        self.primer_alignments_run = True
        for id, seq in enumerate(self.get_sequence_variants()):
            ali_forward = align_primer(seq, primer_forward)
            ali_reverse = align_primer(seq, primer_reverse)
            reverse = id//2
            complementary = id%2
            self.primer_forward_alignments.set_alignment(ali_forward, complementary=complementary, reverse=reverse)
            self.primer_reverse_alignments.set_alignment(ali_reverse, complementary=complementary, reverse=reverse)
    def get_scores(self):
        self.compute_mean_error_rate()
        primer_forward_confidence = self.primer_forward_alignments.get_main_alignment_properties().get("confidence")
        primer_reverse_confidence = self.primer_reverse_alignments.get_main_alignment_properties().get("confidence")
        scores = {"mean_error_rate": self.mean_error_rate,
                  "primer_forward_confidence": primer_forward_confidence,
                  "primer_reverse_confidence": primer_reverse_confidence,
                  "primer_aligned": self.primer_aligned["primer"]}
        return scores
    def print_phred(self):
        """
        Print rescaled phred scores (between 0-9)
        This is useful to print it below a sequence
        -------------------
        score -> error rate |
            0 -> 25%-100%   |
            1 -> 6%-25%     |
            2 -> 1.5%-6%    |
            3 -> 0.4%-1.5%  |
            4 -> 0.1%-0.4%  |
        -------------------
        :return: None
        """
        # for i in range(10):
        #     print(10**((-6*i)/10))
        if self.phred_scores is None:
            print("phred_scores is None")
            return None
        phred_rescaled = np.floor(np.clip(self.phred_scores - 1, 0, 59) / 6).astype(int)  # phred scores rescaled between 0-9
        print("".join([str(phred_rescaled[i]) for i in range(phred_rescaled.size)]))
    def fastq_representation(self):
        """
        Return a fastq string representation of the Sequence instance
        :return: fastq string [str]
        """
        fastq_str = "@" + str(self.sequence_description) + "\n"
        fastq_str = fastq_str + self.sequence + "\n"
        fastq_str = fastq_str + "+\n"
        if self.phred_scores is not None:
            phred_scores = self.phred_scores
        else:
            phred_scores = len(self.sequence) * [0]
        fastq_str = fastq_str + "".join([chr(33 + c) for c in phred_scores]) + "\n"
        return fastq_str
    def __str__(self):
        """
        print(Sequence()) uses this function to format the Sequence() into string representation.
        Useful to print some informations about the Sequence() instance.
        :return: string representation of Sequence() instance [str]
        """
        string = "Sequence() instance\n"
        try:
            string = string + f"Mean error rate: {round(100 * self.mean_error_rate, 2)}%\n"
        except:
            pass
        try:
            primer_forward_alignment = self.primer_forward_alignments.get_formatted_main_alignment(self.sequence)
            properties_forward = self.primer_forward_alignments.get_main_alignment_properties()
            primer_forward_score, primer_forward_start, primer_forward_end = properties_forward.get("score"), properties_forward.get("start"), properties_forward.get("end")
            string = string + f"Primer forward: {primer_forward_alignment[2]} (score={round(primer_forward_score, 2)}, start={primer_forward_start}, end={primer_forward_end})\n"
            string = string + f"                {primer_forward_alignment[1]}\n"
            string = string + f"        {self.sequence[primer_forward_start - 8:primer_forward_start] + primer_forward_alignment[0] + self.sequence[primer_forward_end:primer_forward_end + 8]}\n"

            primer_reverse_alignment = self.primer_reverse_alignments.get_formatted_main_alignment(self.sequence)
            properties_reverse = self.primer_reverse_alignments.get_main_alignment_properties()
            primer_reverse_score, primer_reverse_start, primer_reverse_end = properties_reverse.get("score"), properties_reverse.get("start"), properties_reverse.get("end")
            string = string + f"Primer reverse: {primer_reverse_alignment[2]} (score={round(primer_reverse_score, 2)}, start={primer_reverse_start}, end={primer_reverse_end})\n"
            string = string + f"                {primer_reverse_alignment[1]}\n"
            string = string + f"        {self.sequence[primer_reverse_start - 8:primer_reverse_start] + primer_reverse_alignment[0] + self.sequence[primer_reverse_end:primer_reverse_end + 8]}\n"
        except:
            pass
        return string
class MultiSequence():
    def __init__(self, filepath=None):
        self.sequences = [] # list of Sequence() instances [list of Sequence()]
        self.filepath = filepath # filepath used to load data [str]
        self.consensus = None # consensus sequence [str]
        self.consensus_phred_scores = None # consensus sequence Phred scores [np.ndarray]
        self.run_infos = {"primer_alignments_run": False, # set primer_alignments_run to True when run_primer_alignments() has been run [bool]
                          "ranking_computed": False}
        self.ranking_forward = None
        self.ranking_reverse = None
        if self.filepath is not None:
            self.read_file()
    def read_file(self):
        """
        Read the file in self.filepath, 2 file extension are accepted:
            - .fast: load the data
            - .pickle: load the MultiSequence() previously saved
        /!\ Warning: This function is intended to be used once (for a MultiSequence() instance)
                     It means that multiple .fastq files should be concatenated before being read by this method
        :return: None
        """
        if self.filepath is not None:
            extension = self.filepath[self.filepath.rfind('.'):]
            if extension == ".fastq":
                reads = read_fastq(fastq_filepath=self.filepath)
                for read in reads:
                    self.sequences.append(Sequence(read))
            elif extension == ".pickle":
                self.load(self.filepath)
            else:
                warnings.warn("Error in MultiSequence(): filepath extension not recognised during read_file().")
        else:
            warnings.warn("Error in MultiSequence(): no filepath found during read_file().")
    def write_to_fastq(self, dst=None):
        """
        Write all the sequences in a single .fastq formatted file
        :param dst: filepath of the fastq file destination [str]
        :return: None
        """
        if dst is None: dst = "MultiSequence_saved.fastq"  # default destination (example)
        if dst[-6:] != ".fastq": dst = f"{dst}.fastq"
        with open(dst, mode="w") as f:
            for Seq in self.sequences:
                f.write(Seq.fastq_representation())
    def save(self, dst=None):
        """
        Save the current instance of MultiSequence() to the file dst
        :param dst: filepath of destination [str]
        :return: None
        """
        if dst is None: dst = "MultiSequence_saved.pickle"  # default destination (example)
        if dst[-7:] != ".pickle": dst = f"{dst}.pickle"
        pickle.dump(self.__dict__, open(dst, mode='wb'))
    def load(self, src=None):
        """
        Load a MultiSequence() instance saved
        :param src: filepath of the pickle file source [str]
        :return: None
        """
        if src is None: src = "MultiSequence_saved.pickle"  # default source (example)
        if src[-7:] != ".pickle": src = f"{src}.pickle"
        state = pickle.load(open(src, mode='rb'))
        self.__dict__.update(state)
    def run_primer_alignments(self, primer_forward, primer_reverse):
        """
        Run primer_align on each self.sequences for each Sequence.get_variants() for each primer (forward and reverse)
        This method initializes/sets the primers properties in each Sequence()
        This method prints a progress bar.
        :param primer_forward: 5'->3' forward primer [str]
                               ex: "ATGTCACCACAAACAGAGACTAAAGC"
        :param primer_reverse: 5'->3' reverse primer [str]
                               ex: "TCGCATGTACCTGCAGTAGC"
        :return: None
        """
        time_last = time.time()
        mean_t = 0.0
        for Seq_id in range(len(self.sequences)):
            self.sequences[Seq_id].run_primer_alignments(primer_forward, primer_reverse)
            self.sequences[Seq_id].set_sequence_with_best_primer_alignment()

            delta_t = time.time() - time_last
            time_last = time.time()
            if mean_t == 0: mean_t = delta_t + 1e-10
            mean_t = mean_t * 0.99 + delta_t * 0.01 # "running" time
            print(f"\rrun_primer_alignments progress: {Seq_id+1}/{len(self.sequences)} ({round(1 / mean_t)} reads/sec)", end='')
        print("")
        self.run_infos["primer_alignments_run"] = True
    def compute_sequence_ranking(self):
        sequences_scores = []
        for Seq_id in range(len(self.sequences)):
            scores = self.sequences[Seq_id].get_scores()
            sequences_scores.append((Seq_id, scores))
        ranking_forward = sorted(sequences_scores, key=lambda x: (x[1]["primer_forward_confidence"], x[1]["mean_error_rate"]), reverse=True)
        self.ranking_forward = [r[0] for r in ranking_forward]
        ranking_reverse = sorted(sequences_scores, key=lambda x: (x[1]["primer_reverse_confidence"], x[1]["mean_error_rate"]), reverse=True)
        self.ranking_reverse = [r[0] for r in ranking_reverse]
        self.run_infos["ranking_computed"] == True
    def get_best_sequences(self, n=100, key="forward", skipped=0, truncate_sequences=True):
        """
        Get the `n` best sequences according to the [forward/reverse] primer
        :param n: number of sequences [int]
        :param key: "forward" or "reverse" [str]
        :param skipped: number of best sequences "skipped" before getting n sequences [int]
        :return: n best sequences according to the [forward/reverse] primer
        """
        if self.run_infos["ranking_computed"] == False:
            self.compute_sequence_ranking()
            warnings.warn("[MultiSequence.get_best_sequences] compute_sequence_ranking() has been called to get the best sequences")
        if key == "forward":
            ranking = self.ranking_forward
        elif key == "reverse":
            ranking = self.ranking_reverse
        else:
            warnings.warn(f"[MultiSequence.get_best_sequences] 'key' argument must be equal to either 'forward' or 'reverse'")
        indices = ranking[skipped:skipped+n]
        best_sequences = copy.deepcopy([self.sequences[index] for index in indices])
        if truncate_sequences:
            for Seq_id in range(len(best_sequences)):
                best_sequences[Seq_id].truncate_sequence_using_best_primer_alignment()
        BestSequences = MultiSequence()
        BestSequences.sequences = best_sequences
        return BestSequences
    def apply_threshold_error_rate_max(self, max_error_rate=0.06):
        """
        Apply a threshold (maximum allowed for the mean_error_rate) and remove the Sequence() not in the allowed range
        :param max_error_rate: maximum mean_error_rate allowed
        :return: None
        """
        self.sequences = [seq for seq in self.sequences if seq.get_scores()["mean_error_rate"] <= max_error_rate]
        self.compute_sequence_ranking()
    def apply_threshold_primer_alignment_confidence_min(self, min_primer_alignment_confidence=0.8):
        self.sequences = [seq for seq in self.sequences if seq.get_scores()["primer_alignment_confidence"] >= min_primer_alignment_confidence]
        self.compute_sequence_ranking()
    def print_stats_sequences(self):
        """
        Print informations/stats about the sequences (base composition, sequence lengths)
        :return: None
        """
        if len(self.sequences) == 0:
            warnings.warn("No sequence found in MultiSequence() instance.")
            return None
        seq_lengths = []
        base_A = 0
        base_T = 0
        base_G = 0
        base_C = 0
        for seq in self.sequences:
            seq_lengths.append(len(seq.sequence))
            base_A += seq.sequence.count('A')
            base_T += seq.sequence.count('T')
            base_G += seq.sequence.count('G')
            base_C += seq.sequence.count('C')
        base_total = base_A + base_T + base_G + base_C
        seq_mean_length = np.mean(np.array(seq_lengths))
        seq_std_length = np.std(np.array(seq_lengths))
        seq_median_length = np.median(np.array(seq_lengths))
        print(f"Stats about the sequences (n={len(self.sequences)}):")
        print(f"\tPercentage of A: {round(100 * base_A / base_total, 2)}%")
        print(f"\tPercentage of T: {round(100 * base_T / base_total, 2)}%")
        print(f"\tPercentage of G: {round(100 * base_G / base_total, 2)}%")
        print(f"\tPercentage of C: {round(100 * base_C / base_total, 2)}%")
        print(f"\tMean length:   {round(seq_mean_length, 1)}")
        print(f"\tStd length:    {round(seq_std_length, 1)}")
        print(f"\tMedian length: {round(seq_median_length, 1)}")
        print(f"\tPercentiles length: ", end='')
        for percent in [1, 5, 10, 25, 50, 75, 90, 95, 99]:
            print(f"{round(np.percentile(seq_lengths, percent))} ", end='')
        print(f"\n\t{20 * ' '}1%  5%  10% 25% 50% 75% 90% 95% 99%")
    def print_stats_error_rates(self):
        """
        Print informations/stats about the sequence error rates (sequencing error rates)
        :return: None
        """
        if len(self.sequences) == 0:
            warnings.warn("No sequence found in MultiSequence() instance.")
            return None
        seq_error_rates = []
        for seq in self.sequences:
            mean_error_rate = seq.get_scores()["mean_error_rate"]
            if mean_error_rate is not None: seq_error_rates.append(mean_error_rate)
        if len(seq_error_rates) == 0:
            print("No information about the phred scores")
            return None
        count_error_rate = np.round(100 * np.array(seq_error_rates))
        values, counts = np.unique(count_error_rate, return_counts=True)
        cumulative_percentage = 0.0
        print(f"Stats about the sequencing errors (n={len(self.sequences)}):")
        print(f"\tInfos: number of sequences (proportion|cumulative)")
        for i in range(len(values)):
            percentage = 100 * counts[i] / count_error_rate.size
            cumulative_percentage += percentage
            print(f"\t{round(values[i])}% of errors: {counts[i]} ({int(percentage)}%|{int(cumulative_percentage)}%)")
        error_rate_median = np.median(np.array(seq_error_rates))
        error_rate_mean = np.mean(np.array(seq_error_rates))
        error_rate_std = np.std(np.array(seq_error_rates))
        print(f"\tMedian error rate: {round(100 * error_rate_median, 2)}%")
        print(f"\tMean error rate: {round(100 * error_rate_mean, 2)}%")
        print(f"\tStd error rate: {round(100 * error_rate_std, 2)}%")
    def print_stats(self):
        """
        Print all implemented informations/stats about MultiSequence() instance
        :return: None
        """
        if len(self.sequences) == 0:
            warnings.warn("No sequence found in MultiSequence() instance.")
            return None
        self.print_stats_sequences()
        self.print_stats_error_rates()
    def __str__(self):
        """
        print(MultiSequence()) uses this function to format the MultiSequence() into string representation.
        When printing a MultiSequence() instance: print_stats() is called
        :return: empty string [str]
        """
        print("/!\ Printing a MultiSequence() instance calls the method print_stats().")
        self.print_stats()
        return ""

import unittest
class TestVersion(unittest.TestCase):
    def test_python(self):
        self.assertTrue(sys.version_info >= (3,10)) # Check that Python 3.10 at least is used
    def test_biopython(self):
        import Bio
        self.assertTrue(Bio.__version__ == "1.79", f"biopython version 1.79 needed (version {Bio.__version__} found)")
class TestPrimer(unittest.TestCase):
    def get_Primer(self):
        Primer = PrimerAllAlignments(primer="ATGTCACCACAAACAGAGACTAAAGC")
        Primer.set_alignment(PrimerAlignment(primer_score=2, primer_start=100, primer_end=125), complementary=False, reverse=False)
        Primer.set_alignment(PrimerAlignment(primer_score=13, primer_start=90, primer_end=115), complementary=False, reverse=True)
        Primer.set_alignment(PrimerAlignment(primer_score=20, primer_start=100, primer_end=130), complementary=True, reverse=False)
        Primer.set_alignment(PrimerAlignment(primer_score=40, primer_start=80, primer_end=120), complementary=True, reverse=True)
        return Primer
    def test_flip(self):
        Primer = self.get_Primer()
        Primer.flip_reverse()
        properties = Primer.get_main_alignment_properties()
        self.assertTrue(properties.get("score") == 13)
        self.assertTrue(properties.get("start") == 90)
        self.assertTrue(properties.get("end") == 115)
        Primer.flip_complementary()
        properties = Primer.get_main_alignment_properties()
        self.assertTrue(properties.get("score") == 40)
        self.assertTrue(properties.get("start") == 80)
        self.assertTrue(properties.get("end") == 120)
    def test_resize(self):
        Primer = self.get_Primer()
        Primer.truncate(10,5)
        self.assertEqual(Primer.alignments[0][0].primer_start, 90)
        self.assertEqual(Primer.alignments[0][0].primer_end, 115)
        self.assertEqual(Primer.alignments[0][1].primer_start, 85)
        self.assertEqual(Primer.alignments[0][1].primer_end, 110)
        self.assertEqual(Primer.alignments[1][0].primer_start, 90)
        self.assertEqual(Primer.alignments[1][0].primer_end, 120)
        self.assertEqual(Primer.alignments[1][1].primer_start, 75)
        self.assertEqual(Primer.alignments[1][1].primer_end, 115)
        Primer.extend(10,5)
        self.assertEqual(Primer.alignments[0][0].primer_start, 100)
        self.assertEqual(Primer.alignments[0][0].primer_end, 125)
        self.assertEqual(Primer.alignments[0][1].primer_start, 90)
        self.assertEqual(Primer.alignments[0][1].primer_end, 115)
        self.assertEqual(Primer.alignments[1][0].primer_start, 100)
        self.assertEqual(Primer.alignments[1][0].primer_end, 130)
        self.assertEqual(Primer.alignments[1][1].primer_start, 80)
        self.assertEqual(Primer.alignments[1][1].primer_end, 120)
class TestFunctions(unittest.TestCase):
    def test_damage_sequence(self):
        random.seed(0)
        sequence = 1000*"A"
        damaged = damage_sequence(sequence, mutation_rate=0.0, deletion_rate=0.0, insertion_rate=0.0)
        self.assertTrue(sequence == damaged)
        damaged = damage_sequence(sequence, mutation_rate=1.0, deletion_rate=0.0, insertion_rate=0.0)
        self.assertTrue(damaged.count('A') == 0)
        damaged = damage_sequence(sequence, mutation_rate=0.0, deletion_rate=1.0, insertion_rate=0.0)
        self.assertTrue(len(damaged) == 0)
        damaged = damage_sequence(sequence, mutation_rate=0.0, deletion_rate=0.0, insertion_rate=0.2)
        self.assertTrue(len(damaged) > len(sequence))
def run_tests():
    unittest.main()

if __name__ == '__main__':
    """ Run tests (uncomment this part) """
    # run_tests()

    """ Example of use """
    # concatenate_fastq(src_folder="fastq_pass", dst="rbcL_Qiagen_tomato.fastq") # skip if the fastq are already concatenated in a single file
    MultiSeq = MultiSequence("rbcL_Qiagen_tomato_5000.fastq")
    MultiSeq.run_primer_alignments(primer_forward="ATGTCACCACAAACAGAGACTAAAGC", primer_reverse="TCGCATGTACCTGCAGTAGC") # rbcL primers in this case
    BestSequences = MultiSeq.get_best_sequences(n=25, key="forward", truncate_sequences=True)
    BestSequences.print_stats() # this method can be used before the thresholding to determine the thresholds to use
    BestSequences.save("BestSequences_rbcL_Qiagen_tomato_5000.pickle") # save for later use (without having to rerun run_primer_alignments for example)
    BestSequences = MultiSequence("BestSequences_rbcL_Qiagen_tomato_5000.pickle") # the saved MultiSequence() can be reloaded like this
    BestSequences.write_to_fastq("BestSequences_rbcL_Qiagen_tomato_5000.fastq")

    """ Generate random sequence and damage it (uncomment this part) """
    sequence = create_random_sequence(length=50) # create a random sequence with 50 bases
    sequence_damaged = damage_sequence(sequence, mutation_rate=0.1, deletion_rate=0.1, insertion_rate=0.1) # damage the sequence
    print(f"{sequence}\n{sequence_damaged}")

