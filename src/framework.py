"""
Tester les algorithmes pour le choix d'un algorithme efficace et juste (vitesse et justesse)

Pour augmenter la rapidite, utiliser la parallelisation?

For now, we are doing alignments after concatenating the reads we've got from NP

input: 
    - sequence retrieved from a chosen algorithm
    - Reference sequence
output:
    - alignment score

"""
# from sequence_preprocessing import *
import Bio
from Bio import SeqIO
from sequence_preprocessing import *


def load_fasta(filename):
    for seq_record in SeqIO.parse(filename, "fasta"):
        seq = seq_record.seq
    return seq


class ConsensusSeq:
    def __init__(self, name, filename=None, seq=None):
        self.name = name  # name of the algorithm
        if filename:
            self.filename = filename
            # self.sequence=read_fastq(filename)
            self.sequence = load_fasta(filename)  # output in fastq?
        elif seq:
            self.sequence = seq
        self.length = len(self.sequence)
        # metrics
        self.score = 0  # alignment score
        self.evalue = 0
        self.coverage = 0  # percentage of coverage

    def align_ref(self, ref):
        aligner = Bio.Align.PairwiseAligner()
        alignments = aligner.align(self.sequence, ref)
        self.score = alignments.score
        self.coverage = self.length / len(ref)
        print("Score: " + str(self.score))
        print("Coverage: " + str(self.coverage))


# MultiSeq = MultiSequence("rbcL_MN_tomato.pickle")
ref = load_fasta("rbcL_NCBI_tomato.fasta")
output = ConsensusSeq("test", "rbcL_NCBI_tomato.fasta")
# print(MultiSeq.consensus) #아직 다 안나온듯 함
output.align_ref(ref)
