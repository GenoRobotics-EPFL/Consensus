import unittest
import sys
##sys.path.append("src/sequence_preprocessing_folder")
##from sequence_preprocessing_folder import* 
import MultiSequence 
import PrimerAllAlignments
import PrimerAlignment
import Sequence 
import MultiSequence

import random, sys

class TestVersion(unittest.TestCase):
    def test_python(self):
        self.assertTrue(sys.version_info >= (3,10)) # Check that Python 3.10 at least is used
    def test_biopython(self):
        import Bio
        self.assertTrue(Bio.__version__ == "1.79", f"biopython version 1.79 needed (version {Bio.__version__} found)")
class TestPrimer(unittest.TestCase):

    def get_Primer(self):
        #Primer = spf.PrimerAlignment.PrimerAllAlignments(primer="ATGTCACCACAAACAGAGACTAAAGC")
        Primer = PrimerAllAlignments(primer="ATGTCACCACAAACAGAGACTAAAGC")
        Primer.set_alignment(PrimerAlignment(primer_score=2, primer_start=100, primer_end=125), complementary=False, reverse=False)
        Primer.set_alignment(PrimerAlignment.PrimerAlignment(primer_score=13, primer_start=90, primer_end=115), complementary=False, reverse=True)
        Primer.set_alignment(PrimerAlignment.PrimerAlignment(primer_score=20, primer_start=100, primer_end=130), complementary=True, reverse=False)
        Primer.set_alignment(PrimerAlignment.PrimerAlignment(primer_score=40, primer_start=80, primer_end=120), complementary=True, reverse=True)
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
        damaged = Sequence.damage_sequence(sequence, mutation_rate=0.0, deletion_rate=0.0, insertion_rate=0.0)
        self.assertTrue(sequence == damaged)
        damaged = Sequence.damage_sequence(sequence, mutation_rate=1.0, deletion_rate=0.0, insertion_rate=0.0)
        self.assertTrue(damaged.count('A') == 0)
        damaged = Sequence.damage_sequence(sequence, mutation_rate=0.0, deletion_rate=1.0, insertion_rate=0.0)
        self.assertTrue(len(damaged) == 0)
        damaged = Sequence.damage_sequence(sequence, mutation_rate=0.0, deletion_rate=0.0, insertion_rate=0.2)
        self.assertTrue(len(damaged) > len(sequence))
def run_tests():
    unittest.main()

if __name__ == '__main__':
    """ Run tests (uncomment this part) """
    # run_tests()

    """ Example of use """
    # concatenate_fastq(src_folder="fastq_pass", dst="rbcL_Qiagen_tomato.fastq") # skip if the fastq are already concatenated in a single file
    MultiSeq = MultiSequence.MultiSequence("rbcL_Qiagen_tomato_5000.fastq")
    MultiSeq.run_primer_alignments(primer_forward="ATGTCACCACAAACAGAGACTAAAGC", primer_reverse="TCGCATGTACCTGCAGTAGC") # rbcL primers in this case
    BestSequences = MultiSeq.get_best_sequences(n=25, key="forward", truncate_sequences=True)
    BestSequences.print_stats() # this method can be used before the thresholding to determine the thresholds to use
    BestSequences.save("BestSequences_rbcL_Qiagen_tomato_5000.pickle") # save for later use (without having to rerun run_primer_alignments for example)
    BestSequences = MultiSequence("BestSequences_rbcL_Qiagen_tomato_5000.pickle") # the saved MultiSequence() can be reloaded like this
    BestSequences.write_to_fastq("BestSequences_rbcL_Qiagen_tomato_5000.fastq")

    """ Generate random sequence and damage it (uncomment this part) """
    sequence = Sequence.create_random_sequence(length=50) # create a random sequence with 50 bases
    sequence_damaged = Sequence.damage_sequence(sequence, mutation_rate=0.1, deletion_rate=0.1, insertion_rate=0.1) # damage the sequence
    print(f"{sequence}\n{sequence_damaged}")

