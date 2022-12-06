"""
This version of the file is a TEST version (V1).
It means that some methods are not finished or not optimized.
It also means that everything here is subject to change in later versions.
- Jeremy Goumaz
"""

import gzip, glob, math, random, pickle, time, warnings
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
def align_primer(sequence, primer, print_alignment=False):
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
    aligner.match_score = 3
    aligner.mismatch_score = -2
    aligner.gap_score = -2
    aligner.query_end_gap_score = 0 # must be 0 (to allow gap around the primer for global alignment) [query ~ sequence here]
    aligner.target_end_gap_score = 0 # set to -100 if the primer should be fully integrated in the sequence [target ~ primer here]

    # Alignment
    alignments = aligner.align(sequence, primer)
    alignment = sorted(alignments)[0] # get the best alignment only
    if print_alignment:
        print(f"Score = {round(alignment.score,1)}")
        print(alignment)

    # 'score', 'loc_start', 'loc_end' computations
    score = alignment.score/len(primer)
    loc_start = alignment.aligned[0][0][0]
    loc_end = alignment.aligned[0][-1][1]

    # 'primer_alignment' contains formatted string used to print the alignment results
    aform = alignment._format_pretty()
    aform_linebreak1 = aform.find('\n')
    aform_linebreak2 = aform.find('\n', aform_linebreak1 + 1)
    aform_linebreak3 = aform.find('\n', aform_linebreak2 + 1)
    primer_position_str = aform[aform_linebreak2+1:aform_linebreak3].replace("A", "N").replace("T", "N").replace("G", "N").replace("C", "N")
    pos1 = primer_position_str.find("N")
    pos2 = primer_position_str.rfind("N") + 1
    primer_alignment = []
    primer_alignment.append(aform[0:aform_linebreak1][pos1:pos2])
    primer_alignment.append(aform[aform_linebreak1+1:aform_linebreak2][pos1:pos2])
    primer_alignment.append(aform[aform_linebreak2+1:aform_linebreak3][pos1:pos2])
    return [score, loc_start, loc_end, primer_alignment]
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
    primer = create_random_sequence(primer_length)
    scores = np.empty(shape=(iterations,))
    for it in range(iterations):
        seq = create_random_sequence(seq_length)
        scores[it] = align_primer(seq, primer)[0]
    if print_scores:
        print(f"Random alignments scores (seq_length={seq_length}, primer_length={primer_length}, it={iterations}):")
        for percent in [50, 75, 90, 99, 99.9]:
            print(f"\tBest {round(100-percent,1)}%: {round(np.percentile(scores, percent),2)}")
    return np.percentile(scores, 99.9)

class Sequence():
    def __init__(self, seq=None):
        self.sequence = None  # sequence [str]
        self.phred_scores = None  # Phred scores [np.ndarray]
        self.sequence_description = None  # sequence description from fastq files [str]
        self.mean_error_rate = None  # mean error rate computed from Phred scores [float]
        if seq is not None: self.set_sequence(seq)
        self.primer_score = None # best primer score (max of primer_forward_score and primer_reverse_score) [float]
        self.primer_forward = None  # forward primer sequence [str]
        self.primer_forward_score = None  # alignment score (per base) between forward primer and self.sequence [float]
        self.primer_forward_start = None  # start position of forward primer in self.sequence [int]
        self.primer_forward_end = None  # end position of forward primer in self.sequence [int]
        self.primer_forward_alignment = None  # alignment information formatted (between forward primer and self.sequence) [list of str]
        self.primer_reverse = None  # reverse primer sequence [str]
        self.primer_reverse_score = None  # alignment score (per base) between reverse primer and self.sequence [float]
        self.primer_reverse_start = None  # start position of reverse primer in self.sequence [int]
        self.primer_reverse_end = None  # end position of reverse primer in self.sequence [int]
        self.primer_reverse_alignment = None  # alignment information formatted (between reverse primer and self.sequence) [list of str]
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
                        if id >= 2: # reversed case
                            self.phred_scores =np.flip(self.phred_scores)
                        self.phred_scores = self.phred_scores[loc_start:loc_end]
                if found == False or solvable == False:
                    self.phred_scores = None
            self.sequence = seq
        elif type(seq) is SeqRecord.SeqRecord:
            self.sequence = str(seq.seq)
            self.phred_scores = np.array(seq.letter_annotations['phred_quality'])
            error_proba = 10 ** (-self.phred_scores / 10)
            self.mean_error_rate = error_proba.mean()
            self.sequence_description = str(seq.description)
        else:
            warnings.warn("The sequence type given to Sequence() class is not recognized.")
    def set_phred_scores(self, phred_scores):
        """
        Set Phred scores
        :param phred_scores: phred_scores [np.ndarray | 1D]
        :return: None
        """
        if type(self.phred_scores) != np.ndarray:
            phred_scores = np.array(phred_scores).reshape((-1,))
        self.phred_scores = phred_scores
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
            - sequence
            - complementary sequence
            - reversed sequence
            - reversed complementary sequence
        :return: list of sequences [list of str | len=4]
        """
        sequence_comp = self.get_sequence_complementary()
        return [self.sequence, sequence_comp, self.sequence[::-1], sequence_comp[::-1]]
    def set_forward_primer(self, primer_forward, score=None, start=None, end=None, primer_alignment=None):
        """
        Set the forward primer and its properties
        :param primer_forward: primer sequence [str]
        :param score: alignment score with self.sequence [float]
        :param start: start position in self.sequence [int]
        :param end: end position in self.sequence [int]
        :param primer_alignment: list of 3 formatted strings containing the alignment [list of str | len=3]
            primer_alignment[0]: self.sequence aligned
            primer_alignment[1]: alignments
            primer_alignment[2]: primer sequence aligned
        :return: None
        """
        self.primer_forward = primer_forward
        self.primer_forward_score = score
        self.primer_forward_start = start
        self.primer_forward_end = end
        self.primer_forward_alignment = primer_alignment
        if self.primer_score is None: self.primer_score = -10.0
        if self.primer_forward_score > self.primer_score: self.primer_score = self.primer_forward_score
    def set_reverse_primer(self, primer_reverse, score=None, start=None, end=None, primer_alignment=None):
        """
        Set the reverse primer and its properties
        :param primer_reverse: primer sequence [str]
        :param score: alignment score with self.sequence [float]
        :param start: start position in self.sequence [int]
        :param end: end position in self.sequence [int]
        :param primer_alignment: list of 3 formatted strings containing the alignment [list of str | len=3]
            primer_alignment[0]: self.sequence aligned
            primer_alignment[1]: alignments
            primer_alignment[2]: primer sequence aligned
        :return: None
        """
        self.primer_reverse = primer_reverse
        self.primer_reverse_score = score
        self.primer_reverse_start = start
        self.primer_reverse_end = end
        self.primer_reverse_alignment = primer_alignment
        if self.primer_score is None: self.primer_score = -10.0
        if self.primer_reverse_score > self.primer_score: self.primer_score = self.primer_reverse_score
    def remove_first_bases(self, n=1):
        """
        Remove n nucleotides from the beginning of the sequence
        :param n: number of nucleotides to remove [int]
        :return: None
        """
        if len(self.sequence) != 0: self.sequence = self.sequence[n:]
        if len(self.phred_scores) != 0: self.phred_scores = self.phred_scores[n:]
    def add_first_bases(self, n=1):
        """
        Add n nucleotides (N) in the beginning of the sequence
        :param n: number of nucleotides to add [int]
        :return: None
        """
        for _ in range(n):
            self.sequence = "N" + self.sequence
            self.phred_scores = np.insert(self.phred_scores, 0, 0)
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
        fastq_str = fastq_str + "".join([chr(33 + c) for c in self.phred_scores]) + "\n"
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
            string = string + f"Primer forward: {self.primer_forward_alignment[2]} (score={round(self.primer_forward_score, 2)}, start={self.primer_forward_start}, end={self.primer_forward_end})\n"
            string = string + f"                {self.primer_forward_alignment[1]}\n"
            string = string + f"        {self.sequence[self.primer_forward_start - 8:self.primer_forward_start] + self.primer_forward_alignment[0] + self.sequence[self.primer_forward_end:self.primer_forward_end + 8]}\n"
            string = string + f"Primer reverse: {self.primer_reverse_alignment[2]} (score={round(self.primer_reverse_score, 2)}, start={self.primer_reverse_start}, end={self.primer_reverse_end})\n"
            string = string + f"                {self.primer_reverse_alignment[1]}\n"
            string = string + f"        {self.sequence[self.primer_reverse_start - 8:self.primer_reverse_start] + self.primer_reverse_alignment[0] + self.sequence[self.primer_reverse_end:self.primer_reverse_end + 8]}\n"
        except:
            pass
        return string
class MultiSequence():
    def __init__(self, filepath=None):
        self.sequences = [] # list of Sequence() instances [list of Sequence()]
        self.filepath = filepath # filepath used to load data [str]
        self.primer_forward = None # forward primer sequence [str]
        self.primer_reverse = None # reverse primer sequence [str]
        self.primer_alignments_run = False # set to True when run_primer_alignments() has been run [bool]
        self.consensus = None # consensus sequence [str]
        self.consensus_phred_scores = None # consensus sequence Phred scores [np.ndarray]
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
        :param primer_forward: 5'->3' forward primer [str]
                               ex: "ATGTCACCACAAACAGAGACTAAAGC"
        :param primer_reverse: 5'->3' reverse primer [str]
                               ex: "TCGCATGTACCTGCAGTAGC"
        :return: None
        """
        primer_reverse = primer_reverse.replace("A", "t").replace("T", "a").replace("G", "c").replace("C", "g").upper()[::-1]  # reverse and complementary
        self.primer_forward = primer_forward
        self.primer_reverse = primer_reverse
        self.primer_alignments_run = True
        time_last = time.time()
        mean_t = 0.0
        for Seq_id, Seq in enumerate(self.sequences):
            best_score = 0
            best_ali = [None, None]
            best_seq = None
            for id, seq in enumerate(Seq.get_sequence_variants()):
                ali_forward = align_primer(seq, primer_forward)
                ali_reverse = align_primer(seq, primer_reverse)
                if ali_forward[0] > best_score or ali_reverse[0] > best_score:
                    best_score = max(ali_forward[0], ali_reverse[0])
                    best_ali = [ali_forward, ali_reverse]
                    best_seq = seq
            Seq.set_sequence(best_seq)
            Seq.set_forward_primer(primer_forward, *best_ali[0])
            Seq.set_reverse_primer(primer_reverse, *best_ali[1])

            delta_t = time.time() - time_last
            time_last = time.time()
            if mean_t == 0: mean_t = delta_t
            mean_t = mean_t * 0.99 + delta_t * 0.01 # "running" time
            print(f"\rrun_primer_alignments progress: {Seq_id+1}/{len(self.sequences)} ({round(1 / mean_t)} reads/sec)", end='')
        print("")
    def apply_threshold_error_rate_max(self, max_error_rate=0.06):
        """
        Apply a threshold (maximum allowed for the mean_error_rate) and remove the Sequence() not in the allowed range
        :param max_error_rate: maximum mean_error_rate allowed
        :return: None
        """
        self.sequences = [seq for seq in self.sequences if seq.mean_error_rate <= max_error_rate]
    def apply_threshold_primer_score_min(self, min_primer_score=2.2):
        """
        Apply a threshold (minimum allowed for the primer_score) and remove the Sequence() not in the allowed range
        :param min_primer_score: minimum primer_score allowed
        :return: None
        """
        self.sequences = [seq for seq in self.sequences if seq.primer_score >= min_primer_score]
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
            seq_error_rates.append(seq.mean_error_rate)
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
    def print_stats_primers(self):
        """
        Print informations/stats about the primer alignment scores
        :return: None
        """
        if self.primer_alignments_run == False:
            warnings.warn("run_primer_alignments() has not been run yet.")
            return None
        if len(self.sequences) == 0:
            warnings.warn("No sequence found in MultiSequence() instance.")
            return None
        primer_scores = []
        for seq in self.sequences:
            primer_scores.append(seq.primer_score)
        print(f"Stats about the primers (n={len(self.sequences)}):")
        primer_mean_score = np.mean(np.array(primer_scores))
        primer_std_score = np.std(np.array(primer_scores))
        primer_median_score = np.median(np.array(primer_scores))
        print(f"\tMean score:   {round(primer_mean_score, 1)}")
        print(f"\tStd score:    {round(primer_std_score, 1)}")
        print(f"\tMedian score: {round(primer_median_score, 1)}")
        print(f"\tPercentiles score: ", end='')
        for percent in [1, 5, 10, 25, 50, 75, 90, 95, 99, 99.5, 99.8, 99.9, 99.99]:
            print(f"{str(round(np.percentile(primer_scores, percent),3)).ljust(6)}", end='')
        print(f"\n\t{19 * ' '}1%    5%    10%   25%   50%   75%   90%   95%   99%   99.5% 99.8% 99.9% 99.99%")
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
        self.print_stats_primers()
    def __str__(self):
        """
        print(MultiSequence()) uses this function to format the MultiSequence() into string representation.
        When printing a MultiSequence() instance: print_stats() is called
        :return: empty string [str]
        """
        print("/!\ Printing a MultiSequence() instance calls the method print_stats().")
        self.print_stats()
        return ""



if __name__ == '__main__':
    """ Example of use """
    # concatenate_fastq(src_folder="fastq_pass", dst="rbcL_Qiagen_tomato.fastq") # skip if the fastq are already concatenated in a single file
    MultiSeq = MultiSequence("rbcL_Qiagen_tomato.fastq")
    MultiSeq.run_primer_alignments(primer_forward="ATGTCACCACAAACAGAGACTAAAGC", primer_reverse="TCGCATGTACCTGCAGTAGC") # rbcL primers in this case
    MultiSeq.apply_threshold_primer_score_min(min_primer_score=2.2) # minimum primer alignments score of 2.2
    MultiSeq.apply_threshold_error_rate_max(max_error_rate=0.06) # maximum 6% of errors allowed
    MultiSeq.print_stats() # this method can be used before the thresholding to determine the thresholds to use
    MultiSeq.save("rbcL_Qiagen_tomato.pickle") # save for later use (without having to rerun run_primer_alignments for example)
    MultiSeq = MultiSequence("rbcL_Qiagen_tomato.pickle") # the saved MultiSequence() can be reloaded like this

