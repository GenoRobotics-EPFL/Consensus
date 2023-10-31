"""
This version of the file is a TEST version.
It means that some methods are not finished or not optimized.
It also means that everything here is subject to change in later versions.
- Jeremy Goumaz
"""

__version__ = '2.1.0'

from objects.MultiSequence import MultiSequence
from helpers.fastq import damage_sequence, create_random_sequence
from tests.test import run_tests
from tests.logtime import LogTime, initialize_log_file

if __name__ == '__main__':
    """ Run tests (uncomment this part) """
    # run_tests()

    filename = "assets/input/rbcL_Qiagen_tomato_5000.fastq"
    initialize_log_file(filename)

    """ Example of use """
    # concatenate_fastq(src_folder="fastq_pass", dst="rbcL_Qiagen_tomato.fastq") # skip if the fastq are already concatenated in a single file
    with LogTime("Initialization of MultiSequence"):
        MultiSeq = MultiSequence("assets/input/rbcL_Qiagen_tomato_5000.fastq")

    with LogTime("Running primer alignments"):
        MultiSeq.run_primer_alignments(primer_forward="ATGTCACCACAAACAGAGACTAAAGC", primer_reverse="TCGCATGTACCTGCAGTAGC") # rbcL primers in this case
    
    with LogTime("Getting best sequences"):
        BestSequences = MultiSeq.get_best_sequences(n=25, key="forward", truncate_sequences=True)

    BestSequences.print_stats() # this method can be used before the thresholding to determine the thresholds to use
    BestSequences.save("assets/output/BestSequences_rbcL_Qiagen_tomato_5000.pickle") # save for later use (without having to rerun run_primer_alignments for example)
    BestSequences = MultiSequence("assets/output/BestSequences_rbcL_Qiagen_tomato_5000.pickle") # the saved MultiSequence() can be reloaded like this
    BestSequences.write_to_fastq("assets/output/BestSequences_rbcL_Qiagen_tomato_5000.fastq")

    """ Generate random sequence and damage it (uncomment this part) """
    sequence = create_random_sequence(length=50) # create a random sequence with 50 bases
    sequence_damaged = damage_sequence(sequence, mutation_rate=0.1, deletion_rate=0.1, insertion_rate=0.1) # damage the sequence
    print(f"{sequence}\n{sequence_damaged}")