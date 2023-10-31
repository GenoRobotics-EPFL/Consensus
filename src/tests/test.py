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